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
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <GL/glut.h>
#include <Python.h>

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
#include"PUtils.h"
#include"PM.h"
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

void MainFree(void);
void MainTest(void);
void MainBusyIdle(void);
static void MainInit(void);
void MainReshape(int width, int height);

GLuint obj;


 PyThreadState *_save;

static GLint WinX = 640+cOrthoRightSceneMargin;
static GLint WinY = 480+cOrthoBottomSceneMargin;
static GLint Modifiers = 0;

static char **myArgv;
static int myArgc;

static int FinalInitFlag=1;

int TheWindow;

typedef struct {
  int DirtyFlag;
  int IdleFlag;
  int SwapFlag;
} CMain;

static CMain Main;
int PyMOLReady = false;
int PMGUI = true;

/*========================================================================*/
void MainDirty(void)
{
  CMain *I = &Main;
  I->DirtyFlag=true;
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
  CMain *I = &Main;

  PLock(cLockAPI,&_save);

  /* stay blocked here because Clicks->SetFrame->PParse */

  y=WinY-y;

  glMod = glutGetModifiers();
  Modifiers = ((glMod&GLUT_ACTIVE_SHIFT) ? cOrthoSHIFT : 0) |
	 ((glMod&GLUT_ACTIVE_CTRL) ? cOrthoCTRL : 0) |
	 ((glMod&GLUT_ACTIVE_ALT) ? cOrthoALT : 0);

  if(!OrthoButton(button,state,x,y,Modifiers))
    {
    }
  if(I->SwapFlag)
	 {
		glutSwapBuffers();
		I->SwapFlag=false;
	 }

  PUnlock(cLockAPI,&_save);

}
/*========================================================================*/
static void MainDrag(int x,int y)
{
  CMain *I = &Main;

   PLock(cLockAPI,&_save);

  y=WinY-y;
  if(!OrthoDrag(x,y,Modifiers))
    {
	 }

  if(I->SwapFlag)
	 {
		glutSwapBuffers();
		I->SwapFlag=false;
	 }

   PUnlock(cLockAPI,&_save);

}

/*========================================================================*/
static void MainDraw(void)
{
  CMain *I = &Main;
  int was_dirty;
#ifndef _PYMOL_MODULE
  int a,l;
  char *p;
  char buffer[300];
#endif

  PLock(cLockAPI,&_save);

  was_dirty = I->DirtyFlag;
  if(I->DirtyFlag) I->DirtyFlag=false;
  
  OrthoBusyPrime();
  ExecutiveDrawNow();
  if(FinalInitFlag)
	 {
		Py_BLOCK_THREADS;
		FinalInitFlag=0;

#ifndef _PYMOL_MODULE
		for(a=1;a<myArgc;a++)
		  {
			 l=strlen(myArgv[a]);
			 if(l&&(l<255))
				if(myArgv[a][0]!='-') 
				  {
					 p=strstr(myArgv[a],".py");
					 if(p) {
						sprintf(buffer,"run %s",myArgv[a]);
						PParse(buffer);
					 } else {
						p=strstr(myArgv[a],".pdb");
						if(p) {
						  sprintf(buffer,"load %s",myArgv[a]);
						  PParse(buffer);
						} else {
						  p=strstr(myArgv[a],".mol");
						  if(p) {
							 sprintf(buffer,"load %s",myArgv[a]);
							 PParse(buffer);
						  } else {
							 p=strstr(myArgv[a],".pml");
							 if(p) {
								sprintf(buffer,"@%s",myArgv[a]);
								PParse(buffer);
							 }
						  }
						}
					 }
				  }
		  }
#endif
      
		OrthoRestorePrompt();
		Py_UNBLOCK_THREADS;
	 }
  else if(I->SwapFlag)
    {
      if(PMGUI&&(!was_dirty)&&I->IdleFlag) SceneCopy(0);
      if(PMGUI) glutSwapBuffers();
      I->SwapFlag=false;
    }
  PUnlock(cLockAPI,&_save);
}
/*========================================================================*/
static void MainKey(unsigned char k, int x, int y)
{
  CMain *I = &Main;

  PLock(cLockAPI,&_save);

  switch (k) 
	 {
	 case 27: 
		Py_BLOCK_THREADS;
		PExit(EXIT_SUCCESS);
		break;
	 default:
		OrthoKey(k,x,y);
		break;
	 }

  if(I->SwapFlag)
	 {
		glutSwapBuffers();
		I->SwapFlag=false;
	 }

  PUnlock(cLockAPI,&_save);
  
}

/*========================================================================*/
static void MainSpecial(int k, int x, int y)
{
  CMain *I = &Main;

  PLock(cLockAPI,&_save);

  switch (k)
	 {
	 default:
		break;
	 }
  if(I->SwapFlag)
	 {
		glutSwapBuffers();
		I->SwapFlag=false;
	 }

  PUnlock(cLockAPI,&_save);

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
  /*  GLfloat one[4] = { 1,1,1,1 };*/
  GLfloat low[4] = { 0.20,0.20,0.20,1 };

  CMain *I = &Main;

  I->DirtyFlag=true;
  I->IdleFlag=false;

  if(PMGUI) {
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT,low);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_FALSE);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_NORMALIZE);
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

}

/*========================================================================*/
void MainFree(void)
{
  
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
{ /* should only be called by the master thread */

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
        MainDraw();
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
  int wasDirty = false;

  /* flush command and output queues */

  PLock(cLockAPI,&_save);

  if(ControlIdling()) {
	 SceneIdle(); 
	 I->IdleFlag=false;
  } else {
	 I->IdleFlag=true;
  }

  PFlush(&_save);

  if(I->SwapFlag) {
    if(PMGUI) glutSwapBuffers();
    I->SwapFlag=false;
  }
  if(I->DirtyFlag) {
    wasDirty=true;
    if(PMGUI) 
      glutPostRedisplay();
    else
      MainDraw();
    I->DirtyFlag=false;
  }
  
  if(!wasDirty) {
    if(I->IdleFlag) { /* select to avoid racing the CPU */
      
      PUnlock(cLockAPI,&_save);
      PSleep(20000);
      PLock(cLockAPI,&_save);

      if(I->SwapFlag) {
        if(PMGUI) glutSwapBuffers();
        I->SwapFlag=false;
      }
      if(I->DirtyFlag) {
        if(PMGUI) 
          glutPostRedisplay();
        else
          MainDraw();
        I->DirtyFlag=false;
      }
    }
  }
  
  PUnlock(cLockAPI,&_save);
}

  CMain *I = &Main;

/*========================================================================*/

#ifndef _PYMOL_MODULE
int main(int argc, char *argv[])
{
  int stereo;
  myArgc=argc;
  myArgv=argv;
  
#else
void was_main(void) 
{
  int argc = 1;
  char *argv[2],argvv[1024] = "pymol";
  argv[0]=argvv;
  argv[1]=NULL;

#ifdef _DRI_WORKAROUND
  dlopen("libGL.so.1",RTLD_LAZY|RTLD_GLOBAL);
#endif

#endif  


  myArgc=argc;
  myArgv=argv;

  if(myArgc>1)
    if(strcmp(myArgv[1],"-c")==0)
      PMGUI=false;
  
  if(PMGUI) {
    
    glutInit(&argc, argv);
    
#ifdef _PYMOL_STEREO
    glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE | GLUT_STEREO );
#else
    glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE );      
#endif

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
  }

  Py_UNBLOCK_THREADS;
  
  PyMOLReady = true;

  if(PMGUI) {
    printf(" GL based graphics front end:\n");
    printf("  GL_VENDOR: %s\n",(char*)glGetString(GL_VENDOR));
    printf("  GL_RENDERER: %s\n",(char*)glGetString(GL_RENDERER));
    printf("  GL_VERSION: %s\n",(char*)glGetString(GL_VERSION));
    printf("  GL_EXTENSIONS: %s\n",(char*)glGetString(GL_EXTENSIONS));
    glutMainLoop();
  } else {
    printf(" No graphics front end.\n");
    MainReshape(WinX,WinY);
    MainDraw(); /* for command line processing */
    while(1) MainBusyIdle();
  }

#ifndef _PYMOL_MODULE
  return 0;
#endif

}






