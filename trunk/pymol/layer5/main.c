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

void MainFree(void);
void MainTest(void);
void MainBusyIdle(void);
static void MainInit(void);
void MainReshape(int width, int height);
static void MainDrawLocked(void);

GLuint obj;

static GLint WinX = 640;
static GLint WinY = 480;
static GLint Modifiers = 0;

static char **myArgv,*myArgvv[2],myArgvvv[1024];
static int myArgc;

static int FinalInitFlag=1;

int TheWindow;

typedef struct {
  int DirtyFlag;
  int IdleMode;
  int SwapFlag;
  float IdleTime;
  int IdleCount;
  int ReshapeFlag;
} CMain;

static CMain Main;
int PyMOLReady = false;
int PyMOLTerminating = false;
int PMGUI = true;
int StereoCapable=false;
static int InternalGUI = true;
static int InternalFeedback = true;
int ShowSplash=true;

void launch(void);

void MainOnExit(void);

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
    exit(0);
  }
}
/*========================================================================*/
int MainSavingUnderWhileIdle(void)
{
  CMain *I = &Main;
  return(I->IdleMode==2);
}
/*========================================================================*/
void MainDirty(void)
{
  CMain *I = &Main;
  PRINTFD(FB_Main)
    " MainDirty: called.\n"
    ENDFD;
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
void MainTest(void)
{
}
/*========================================================================*/
static void MainButton(int button,int state,int x,int y)
{
  static int glMod;  
  /*  CMain *I = &Main;*/

  PLockAPIAsGlut();

  /* stay blocked here because Clicks->SetFrame->PParse */

  y=WinY-y;

  glMod = p_glutGetModifiers();
  Modifiers = ((glMod&P_GLUT_ACTIVE_SHIFT) ? cOrthoSHIFT : 0) |
	 ((glMod&P_GLUT_ACTIVE_CTRL) ? cOrthoCTRL : 0) |
	 ((glMod&P_GLUT_ACTIVE_ALT) ? cOrthoALT : 0);

  if(!OrthoButton(button,state,x,y,Modifiers))
    {
    }
  PUnlockAPIAsGlut();

}
/*========================================================================*/
static void MainDrag(int x,int y)
{
  /*  CMain *I = &Main;*/
  
  PLockAPIAsGlut();
  
  y=WinY-y;
  if(!OrthoDrag(x,y,Modifiers))
    {
	 }
  
  PUnlockAPIAsGlut();

}
/*========================================================================*/
static void MainDrawLocked(void)
{
  CMain *I = &Main;

  if(I->DirtyFlag) {
    I->DirtyFlag=false;
  }
  
  OrthoBusyPrime();
  ExecutiveDrawNow();
  if(I->SwapFlag)
    {
      if(!SettingGet(cSetting_suspend_updates))
        if(PMGUI) p_glutSwapBuffers();
      I->SwapFlag=false;
    }
}
/*========================================================================*/
static void MainDraw(void)
{
  PLockAPIAsGlut();
  MainDrawLocked();
  PUnlockAPIAsGlut();
}
/*========================================================================*/
static void MainKey(unsigned char k, int x, int y)
{
  /*  CMain *I = &Main;*/
  int glMod;

  PLockAPIAsGlut();

  glMod = p_glutGetModifiers();
  Modifiers = ((glMod&P_GLUT_ACTIVE_SHIFT) ? cOrthoSHIFT : 0) |
	 ((glMod&P_GLUT_ACTIVE_CTRL) ? cOrthoCTRL : 0) |
	 ((glMod&P_GLUT_ACTIVE_ALT) ? cOrthoALT : 0);

  OrthoKey(k,x,y,Modifiers);

  PUnlockAPIAsGlut();
  
}

/*========================================================================*/
static void MainSpecial(int k, int x, int y)
{
  char buffer[255];
  int grabbed = false;
  PLockAPIAsGlut();
  switch(k) {
    case P_GLUT_KEY_UP:
    case P_GLUT_KEY_DOWN:
      grabbed=1;
      OrthoSpecial(k,x,y);
      break;
    case P_GLUT_KEY_LEFT:
    case P_GLUT_KEY_RIGHT:      
      if(OrthoArrowsGrabbed()) {
        grabbed=1;
        OrthoSpecial(k,x,y);
      }
      break;
  }
  if(!grabbed) {
    sprintf(buffer,"_special %d,%d,%d ",k,x,y);
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
  

  OrthoReshape(width,height);
}
/*========================================================================*/
void MainDoReshape(int width, int height) /* called internally */
{
  int h,w;
  int internal_feedback;
  if(width<0) {
    BlockGetSize(SceneGetBlock(),&width,&h);
    if(SettingGet(cSetting_internal_gui))
      width+=SettingGet(cSetting_internal_gui_width);
  }
  if(height<0) { 
    BlockGetSize(SceneGetBlock(),&w,&height);
    internal_feedback = SettingGet(cSetting_internal_feedback);
    if(internal_feedback)
      height+=(internal_feedback-1)*cOrthoLineHeight+cOrthoBottomSceneMargin;
  }
  if(PMGUI) {
    p_glutReshapeWindow(width,height);
    glViewport(0, 0, (GLint) width, (GLint) height);
  }
  OrthoReshape(width,height);
  if(SettingGet(cSetting_full_screen))
    p_glutFullScreen();

}
/*========================================================================*/
static void MainInit(void) 
{
  /*  GLfloat one[4] = { 1,1,1,1 }; 
      GLfloat low[4] = { 0.20,0.20,0.20,1 };*/

  CMain *I = &Main;

  I->DirtyFlag=true;
  I->IdleMode=2;
  I->IdleTime=UtilGetSeconds();
  I->IdleCount = 0;
  I->ReshapeFlag=false;
  if(PMGUI) {

    /* get us into a well defined GL state */
    
    glDisable(GL_LIGHTING);
    glDisable(GL_FOG);
    glDisable(GL_NORMALIZE);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_LINE_SMOOTH);
    glDisable(GL_DITHER);
    glDisable(GL_BLEND);
    /*    glDisable(GL_ALPHA_TEST);
    glDisable(GL_CULL_FACE);
    glDisable(GL_POINT_SMOOTH);*/
    

/*    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT,low);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_FALSE);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_NORMALIZE);*/
  }

  FeedbackInit();
  UtilInit();
  SettingInitGlobal();  
  SettingSet(cSetting_internal_gui,InternalGUI);
  SettingSet(cSetting_internal_feedback,InternalFeedback);
  SphereInit();
  ColorInit();
  OrthoInit(ShowSplash);
  SelectorInit();
  MovieInit();
  SceneInit();
  SculptCacheInit();
  ExecutiveInit();
  RepMeshInit();
  EditorInit();  

}

/*========================================================================*/
void MainFree(void)
{
  PyMOLTerminating=true;

  EditorFree();
  ExecutiveFree();
  SculptCacheFree();
  SceneFree();
  MovieFree();
  SelectorFree();
  OrthoFree();
  SettingFreeGlobal();
  ColorFree();
  SphereDone();
  PFree();
  FeedbackFree();

  MemoryDebugDump();
  printf(" PyMOL: normal program termination.\n");
  
#ifdef WIN32
  ExitProcess(0); /* Py_Exit hangs on WIN32 for some reason */
#endif

}
/*========================================================================*/
void MainRefreshNow(void) 
{ /* should only be called by the master thread, with a locked API */

  CMain *I = &Main;
  if(I->SwapFlag)
    {
      if(PMGUI) p_glutSwapBuffers();
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

void MainBusyIdle(void) 
{
  /* This is one of the few places in the program where we can be sure 
	* that we have the "glut" thread...glut doesn't seem to be completely
	* thread safe or rather thread consistent
   */

  CMain *I = &Main;
  /* flush command and output queues */
  
  PRINTFD(FB_Main)
    " MainBusyIdle: entered, IdleMode %d, DirtyFlag %d, SwapFlag %d\n",
    I->IdleMode,I->DirtyFlag,I->SwapFlag
    ENDFD;
  PLockAPIAsGlut();

  if(ControlIdling()) {
    ExecutiveSculptIterateAll();
	 SceneIdle(); 
	 I->IdleMode=0;
  } else {
    if(!I->IdleMode) {
      I->IdleTime=UtilGetSeconds();
      I->IdleMode=1;
    }
  }

  PFlush();

  if(I->SwapFlag) {
    if(PMGUI) p_glutSwapBuffers();
    I->SwapFlag=false;
  }
  if(I->DirtyFlag) {
    if(PMGUI) 
      p_glutPostRedisplay();
    else
      MainDrawLocked();
    I->DirtyFlag=false;
  }

  if(I->IdleMode) { /* avoid racing the CPU */
    if(I->IdleMode==1) {
      if(UtilGetSeconds()-I->IdleTime>SettingGet(cSetting_idle_delay)) { 
        I->IdleMode=2;
        if(PMGUI)
          if(SettingGet(cSetting_cache_display))
             p_glutPostRedisplay(); /* trigger caching of the current scene */
      }
    }
    if(I->IdleMode==1)
      PSleep(SettingGet(cSetting_fast_idle)); /* fast idle - more responsive */
    else
      PSleep(SettingGet(cSetting_slow_idle)); /* slow idle - save CPU cycles */
  } else {
    PSleep(SettingGet(cSetting_no_idle)); /* give Tcl/Tk a chance to run */
  }
  PUnlockAPIAsGlut();

  if(FinalInitFlag)
	 {
      FinalInitFlag=FinalInitFlag+1;
      if(FinalInitFlag>=10) {
        FinalInitFlag=0;
        PBlock();
#ifndef _PYMOL_WX_GLUT
        PRunString("launch_gui()");
#endif
        PRunString("adapt_to_hardware()");
        PRunString("exec_deferred()");
        PUnblock();
      }
    }
  if(I->ReshapeFlag) {
    MainDoReshape(WinX,WinY);
    I->ReshapeFlag=false;
  }
  if(!PMGUI) {
    if(!OrthoCommandWaiting()) {
      I->IdleCount++;
      if(I->IdleCount==10) {
        PLockAPIAsGlut();
        PParse("_quit");
        PFlush();
        PUnlockAPIAsGlut();
      }
    }
      
  }
  PRINTFD(FB_Main)
    " MainBusyIdle: leaving... IdleMode %d, DirtyFlag %d, SwapFlag %d\n",
    I->IdleMode,I->DirtyFlag,I->SwapFlag
    ENDFD;

}

/*========================================================================*/

void launch(void)
{
  if(InternalGUI)
    WinX+=cOrthoRightSceneMargin;
  if(InternalFeedback)
    WinY+= (InternalFeedback-1)*cOrthoLineHeight + cOrthoBottomSceneMargin;

  if(PMGUI) {
    
    atexit(MainOnExit); /* register callback to help prevent crashes
                                 when GLUT spontaneously kills us */
  
    p_glutInit(&myArgc, myArgv);

    p_glutInitDisplayMode(P_GLUT_RGBA | P_GLUT_DEPTH | P_GLUT_DOUBLE | P_GLUT_STEREO );
    if(!p_glutGet(P_GLUT_DISPLAY_MODE_POSSIBLE)) {
      p_glutInitDisplayMode(P_GLUT_RGBA | P_GLUT_DEPTH | P_GLUT_DOUBLE );            
      StereoCapable = 0;
    } else {
      StereoCapable = 1;
    }

    p_glutInitWindowPosition(0, 175);
    p_glutInitWindowSize(WinX, WinY);
    
    TheWindow = p_glutCreateWindow("PyMOL Viewer");

  } 

  MainInit();

  PInit();

  if(PMGUI) {
    p_glutDisplayFunc(         MainDraw );
    p_glutReshapeFunc(         MainReshape );
    p_glutKeyboardFunc(        MainKey );
    p_glutMouseFunc(           MainButton );
    p_glutMotionFunc(          MainDrag );
    /*  p_glutPassiveMotionFunc(   MainMove );*/
    p_glutSpecialFunc(         MainSpecial );
    p_glutIdleFunc(         MainBusyIdle );

    p_glutPostRedisplay();
  }

  PUnblock();

  PyMOLReady = true;

  if(PMGUI) {
    SceneSetCardInfo((char*)glGetString(GL_VENDOR),
                     (char*)glGetString(GL_RENDERER),
                     (char*)glGetString(GL_VERSION));
    
    printf(" OpenGL based graphics front end:\n");
    printf("  GL_VENDOR: %s\n",(char*)glGetString(GL_VENDOR));
    printf("  GL_RENDERER: %s\n",(char*)glGetString(GL_RENDERER));
    printf("  GL_VERSION: %s\n",(char*)glGetString(GL_VERSION));
    /*        printf("  GL_EXTENSIONS: %s\n",(char*)glGetString(GL_EXTENSIONS));*/
    if(StereoCapable) {
      printf(" Hardware stereo capability detected.\n");
    } 
    p_glutMainLoop();
    PBlock(); /* if we've gotten here, then we're heading back to Python... */
  } else {
    SceneSetCardInfo("none","ray trace only","none");
    printf(" Command mode. No graphics front end.\n");
    MainReshape(WinX,WinY);
    MainDraw(); /* for command line processing */
    while(1) {
      MainBusyIdle();
      MainDraw();
    }
  }
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

  PGetOptions(&PMGUI,&InternalGUI,&ShowSplash,&InternalFeedback);
  launch();

  return 0;

}






