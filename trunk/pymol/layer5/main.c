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
#include <Python.h>

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
#include"Object.h"
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

void MainFree(void);
void MainTest(void);
void MainBusyIdle(void);
static void MainInit(void);
void MainReshape(int width, int height);
static void MainDrawLocked(void);

GLuint obj;



static GLint WinX = 640+cOrthoRightSceneMargin;
static GLint WinY = 480+cOrthoBottomSceneMargin;
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
} CMain;

static CMain Main;
int PyMOLReady = false;
int PMGUI = true;
int StereoCapable=false;
int InternalGUI=true;

void launch(void);

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

  glMod = glutGetModifiers();
  Modifiers = ((glMod&GLUT_ACTIVE_SHIFT) ? cOrthoSHIFT : 0) |
	 ((glMod&GLUT_ACTIVE_CTRL) ? cOrthoCTRL : 0) |
	 ((glMod&GLUT_ACTIVE_ALT) ? cOrthoALT : 0);

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
      if(PMGUI) glutSwapBuffers();
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

  glMod = glutGetModifiers();
  Modifiers = ((glMod&GLUT_ACTIVE_SHIFT) ? cOrthoSHIFT : 0) |
	 ((glMod&GLUT_ACTIVE_CTRL) ? cOrthoCTRL : 0) |
	 ((glMod&GLUT_ACTIVE_ALT) ? cOrthoALT : 0);

  switch (k) 
	 {
	 case 27: 
      PLockAPIAsGlut();
      PParse("_quit");
      PFlush();
      PUnlockAPIAsGlut();
		break;
	 default:
		OrthoKey(k,x,y,Modifiers);
		break;
	 }

  PUnlockAPIAsGlut();
  
}

/*========================================================================*/
static void MainSpecial(int k, int x, int y)
{
  char buffer[255];
  int grabbed = false;
  PLockAPIAsGlut();
  switch(k) {
    case GLUT_KEY_UP:
    case GLUT_KEY_DOWN:
      grabbed=1;
      OrthoSpecial(k,x,y);
      break;
    case GLUT_KEY_LEFT:
    case GLUT_KEY_RIGHT:      
      if(OrthoArrowsGrabbed()) {
        grabbed=1;
        OrthoSpecial(k,x,y);
      }
      break;
  }
  if(!grabbed) {
    sprintf(buffer,"_special %d,%d,%d ",k,x,y);
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

  if(PMGUI) {
    glutReshapeWindow(width,height);
    glViewport(0, 0, (GLint) width, (GLint) height);
  }

  OrthoReshape(width,height);
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

  UtilInit();
  SettingInit();  
  SphereInit();
  ColorInit();
  OrthoInit();
  SelectorInit();
  MovieInit();
  SceneInit();
  ExecutiveInit();
  RepMeshInit();
  EditorInit();  

}

/*========================================================================*/
void MainFree(void)
{
  
  EditorFree();
  ExecutiveFree();
  SceneFree();
  MovieFree();
  SelectorFree();
  OrthoFree();
  SettingFree();
  ColorFree();
  SphereDone();
  PFree();
  
  MemoryDebugDump();
}
/*========================================================================*/
void MainRefreshNow(void) 
{ /* should only be called by the master thread, with a locked API */

  CMain *I = &Main;
  if(I->SwapFlag)
    {
      if(PMGUI) glutSwapBuffers();
      I->SwapFlag=false;
    }
  if(I->DirtyFlag)
    {
      if(PMGUI) 
        glutPostRedisplay();
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

  PLockAPIAsGlut();

  if(ControlIdling()) {
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
    if(PMGUI) glutSwapBuffers();
    I->SwapFlag=false;
  }
  if(I->DirtyFlag) {
    if(PMGUI) 
      glutPostRedisplay();
    else
      MainDrawLocked();
    I->DirtyFlag=false;
  }
  
  if(I->IdleMode) { /* avoid racing the CPU */
    if(I->IdleMode==1) {
      if(UtilGetSeconds()-I->IdleTime>SettingGet(cSetting_idle_delay)) { 
        I->IdleMode=2;
        if(PMGUI) 
          glutPostRedisplay(); /* trigger caching of the current scene */
      }
    }
    PUnlockAPIAsGlut();
    if(I->IdleMode==1)
      PSleep(SettingGet(cSetting_fast_idle)); /* fast idle - more responsive */
    else
      PSleep(SettingGet(cSetting_slow_idle)); /* slow idle - save CPU cycles */
    PLockAPIAsGlut();
  } else {
    PUnlockAPIAsGlut();
    PSleep(SettingGet(cSetting_no_idle)); /* give Tcl/Tk a chance to run */
    PLockAPIAsGlut();
  }
  PUnlockAPIAsGlut();

  if(FinalInitFlag)
	 {
      FinalInitFlag=FinalInitFlag+1;
      if(FinalInitFlag>=10) {
        FinalInitFlag=0;
        PBlock();
        PRunString("launch_gui()");
        PRunString("exec_deferred()");
        PUnblock();
      }
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
}

/*========================================================================*/

void launch(void)
{
  if(PMGUI) {
    
    glutInit(&myArgc, myArgv);

    if(StereoCapable) {
      glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE | GLUT_STEREO );
    } else {
      glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE );      
    }

    glutInitWindowPosition(0, 175);
    glutInitWindowSize(WinX, WinY);
    
    TheWindow = glutCreateWindow("PyMol Viewer");
  }

  MainInit();

  PInit();

  if(PMGUI) {
    glutDisplayFunc(         MainDraw );
    glutReshapeFunc(         MainReshape );
    glutKeyboardFunc(        MainKey );
    glutMouseFunc(           MainButton );
    glutMotionFunc(          MainDrag );
    /*  glutPassiveMotionFunc(   MainMove );*/
    glutSpecialFunc(         MainSpecial );
    glutIdleFunc(         MainBusyIdle );

    glutPostRedisplay();

    if(StereoCapable) SettingSet(cSetting_line_smooth,1.0);
  }

  PUnblock();

  PyMOLReady = true;

  if(PMGUI) {
    printf(" GL based graphics front end:\n");
    printf("  GL_VENDOR: %s\n",(char*)glGetString(GL_VENDOR));
    printf("  GL_RENDERER: %s\n",(char*)glGetString(GL_RENDERER));
    printf("  GL_VERSION: %s\n",(char*)glGetString(GL_VERSION));
    /*    printf("  GL_EXTENSIONS: %s\n",(char*)glGetString(GL_EXTENSIONS));*/
    glutMainLoop();
  } else {
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

  PGetOptions(&PMGUI,&InternalGUI,&StereoCapable);

  launch();

  return 0;

}






