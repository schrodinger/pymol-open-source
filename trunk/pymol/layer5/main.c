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
#include"Isosurf.h"
#include"Tetsurf.h"
#include"PConv.h"
#include"VFont.h"

void MainFree(void);
void MainTest(void);
void MainBusyIdle(void);
static void MainInit(void);
void MainReshape(int width, int height);
static void MainDrawLocked(void);

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
  float IdleTime;
  int IdleCount;
  int ReshapeFlag;
} CMain;

static CMain Main;
int PyMOLReady = false;
int PyMOLTerminating = false;
int PMGUI = true;
int StereoCapable=false;
int Security = true;
int ForceStereo = 0; /* 1 = force stereo (if possible); -1 = force mono; 0 = autodetect */
int GameMode = false;
int BlueLine = false;

static int InternalGUI = true;
static int InternalFeedback = true;
int ShowSplash=true;

void launch(void);

void MainOnExit(void);

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
      glVertex3f(0.0f, window_height - 0.5f, 0.0f);
      glVertex3f(window_width, window_height - 0.5f, 0.0f);
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

void MainRunString(char *str);
PyObject *MainRunStringBlocked(char *str);
void MainRunCommand(char *str1);

/*========================================================================*/
void MainRunCommand(char *str1)
{
  if(str1[0]!='_') { /* suppress internal call-backs */
    if(strncmp(str1,"cmd._",5)) {
      OrthoAddOutput("PyMOL>");
      OrthoAddOutput(str1);
      OrthoNewLine(NULL,true);
      if(WordMatch(str1,"quit",true)==0) /* don't log quit */
        PLog(str1,cPLog_pml);
    }
    PParse(str1);
  } else if(str1[1]==' ') { /* "_ command" suppresses echoing of command, but it is still logged */
    if(WordMatch(str1+2,"quit",true)==0) /* don't log quit */
      PLog(str1+2,cPLog_pml);
    PParse(str1+2);    
  } else {
    PParse(str1);
  }
}
/*========================================================================*/
void MainRunString(char *str)
{
  PBlock();
  PRunString(str);
  PUnblock();
}
/*========================================================================*/
PyObject *MainRunStringBlocked(char *str)
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
  PLockAPIAsGlut();
  MainDrawLocked();
  PUnlockAPIAsGlut();
}
/*========================================================================*/
static void MainKey(unsigned char k, int x, int y)
{
  /*  CMain *I = &Main;*/
  int glMod;

  PRINTFD(FB_Main)
    " MainKey: %d %d %d\n",k,x,y
    ENDFD;
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

PyObject *MainAsPyList(void) 
{
  PyObject *result=NULL;
  int width,height;
  result = PyList_New(2);
  BlockGetSize(SceneGetBlock(),&width,&height);
  PyList_SetItem(result,0,PyInt_FromLong(width));
  PyList_SetItem(result,1,PyInt_FromLong(height));
  return(PConvAutoNone(result));
}

int MainFromPyList(PyObject *list)
{
  int ok=true;
  int win_x,win_y;
  int size;
  OrthoLineType buffer;
  if(ok) ok = (list!=NULL);
  if(ok) ok = PyList_Check(list);
  if(ok) size=PyList_Size(list);
  if(ok&&(size>=2)) {
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
    glDisable(GL_POLYGON_SMOOTH);
    glDisable(GL_DITHER);
    glDisable(GL_BLEND);
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
  VFontInit();
  ExecutiveInit();
  IsosurfInit();
  TetsurfInit();
  EditorInit();  

}

/*========================================================================*/
void MainFree(void)
{
  PyMOLTerminating=true;

  EditorFree();
  ExecutiveFree();
  VFontFree();
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

  if(ShowSplash) {
    MemoryDebugDump();
    printf(" PyMOL: normal program termination.\n");
  }
  
#ifdef WIN32
  TerminateProcess(GetCurrentProcess(),0); /* only way to avoid a crash */
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

void MainBusyIdle(void) 
{
  /* This is one of the few places in the program where we can be sure 
	* that we have the "glut" thread...glut doesn't seem to be completely
	* thread safe or rather thread consistent
   */

  CMain *I = &Main;
  /* flush command and output queues */
  
  /*  PRINTFD(FB_Main)
    " MainBusyIdle: entered, IdleMode %d, DirtyFlag %d, SwapFlag %d\n",
    I->IdleMode,I->DirtyFlag,I->SwapFlag
    ENDFD;*/
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
    if(PMGUI) {
      DrawBlueLine();
      p_glutSwapBuffers();
    }
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
        /* restore working directory if asked to */
        PRunString("if os.environ.has_key('PYMOL_WD'): os.chdir(os.environ['PYMOL_WD'])");

        PRunString("launch_gui()");
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
  /*  PRINTFD(FB_Main)
    " MainBusyIdle: leaving... IdleMode %d, DirtyFlag %d, SwapFlag %d\n",
    I->IdleMode,I->DirtyFlag,I->SwapFlag
    ENDFD;*/

}

/*========================================================================*/

void launch(void)
{
  if(InternalGUI&&(!GameMode))
    WinX+=cOrthoRightSceneMargin;
  if(InternalFeedback&&(!GameMode))
    WinY+= (InternalFeedback-1)*cOrthoLineHeight + cOrthoBottomSceneMargin;

  if(PMGUI) {
    #ifndef _PYMOL_OSX
    atexit(MainOnExit); /* register callback to help prevent crashes
                                 when GLUT spontaneously kills us */
    #endif

    p_glutInit(&myArgc, myArgv);

    switch(ForceStereo) {

    case -1: /* force mono */
      p_glutInitDisplayMode(P_GLUT_RGBA | P_GLUT_DEPTH | P_GLUT_DOUBLE );
      StereoCapable = 0;
      break;

    case 0: /* default/autodetect (stereo on win/unix; mono on macs) */
#ifndef _PYMOL_OSX
      p_glutInitDisplayMode(P_GLUT_RGBA | P_GLUT_DEPTH | P_GLUT_DOUBLE | P_GLUT_STEREO );
      if(!p_glutGet(P_GLUT_DISPLAY_MODE_POSSIBLE)) {
#endif
        p_glutInitDisplayMode(P_GLUT_RGBA | P_GLUT_DEPTH | P_GLUT_DOUBLE );            
        StereoCapable = 0;
#ifndef _PYMOL_OSX
      } else {
        StereoCapable = 1;
      }
#endif
      break;

    case 1: /* force stereo (if possible) */
      p_glutInitDisplayMode(P_GLUT_RGBA | P_GLUT_DEPTH | P_GLUT_DOUBLE | P_GLUT_STEREO );
      if(!p_glutGet(P_GLUT_DISPLAY_MODE_POSSIBLE)) {
        p_glutInitDisplayMode(P_GLUT_RGBA | P_GLUT_DEPTH | P_GLUT_DOUBLE );            
        StereoCapable = 0;
      } else {
        StereoCapable = 1;
      }      
      break;
    }

    if(!GameMode) {
      /*
	#ifdef _PYMOL_OSX
	p_glutInitWindowPosition(0, 200);
	#else
	p_glutInitWindowPosition(0, 175);
	#endif
      */
      p_glutInitWindowPosition(WinPX,WinPY);
      p_glutInitWindowSize(WinX, WinY);

      TheWindow = p_glutCreateWindow("PyMOL Viewer");
    } else {
      char str[255];
      sprintf(str,"%dx%d:32@120",WinX,WinY);
      glutGameModeString(str);
      glutEnterGameMode(); 
    }
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
    p_glutIdleFunc(            MainBusyIdle );

    p_glutPostRedisplay();
  }

  PUnblock();

  if(PMGUI) {
    SceneSetCardInfo((char*)glGetString(GL_VENDOR),
                     (char*)glGetString(GL_RENDERER),
                     (char*)glGetString(GL_VERSION));
    if(ShowSplash) {
      
      printf(" OpenGL based graphics front end:\n");
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
    p_glutMainLoop();
    PBlock(); /* if we've gotten here, then we're heading back to Python... */
  } else {
    SceneSetCardInfo("none","ray trace only","none");
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

  PGetOptions(&PMGUI,&InternalGUI,&ShowSplash,
	      &InternalFeedback,&Security,&GameMode,
	      &ForceStereo,&WinX,&WinY,&BlueLine,
	      &WinPX,&WinPY);
  launch();

  return 0;

}






