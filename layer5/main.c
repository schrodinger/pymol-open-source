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

GLuint obj;

int TheWindow;

 PyThreadState *_save;

static GLint WinX = 640+cOrthoRightSceneMargin;
static GLint WinY = 480+cOrthoBottomSceneMargin;
static GLint Modifiers = 0;

static char **myArgv;
static int myArgc;

static int FinalInitFlag=1;

typedef struct {
  int DirtyFlag;
  int IdleFlag;
  int SwapFlag;
} CMain;

static CMain Main;
int PyMOLReady = false;

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
#ifndef _PYMOL_MODULE
  int a,l;
  char *p;
  char buffer[300];
#endif

  PLock(cLockAPI,&_save);

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

  if(I->SwapFlag)
	 {
		glutSwapBuffers();
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
void MainReshape(int width, int height)
{
  float h;

  WinX = width;
  WinY = height;

  h = ((float)height)/width;
  glViewport(0, 0, (GLint) width, (GLint) height);
  /*  glMatrixMode(GL_PROJECTION);*/
  /* glLoadIdentity(); */
  /*  glFrustum(-1.0, 1.0, -h, h, 5.0, 60.0);*/
  /*  glMatrixMode(GL_MODELVIEW);*/
  /*  glLoadIdentity();*/
  /*  glTranslatef(0.0, 0.0, -40.0); */
  glutReshapeWindow(width, height);
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

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT,low);
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_FALSE);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_NORMALIZE);

  UtilInit();
  SettingInit();  
  SphereInit();
  ColorInit();
  OrthoInit();
  /*  MenuInit();*/
  SelectorInit();
  MovieInit();
  SceneInit();
  ExecutiveInit();

  RepMeshInit();
  /*  CPKInit();*/

  /*  MenuAdd(0,"File",NULL);
  MenuAdd(0,"Quit",MainExit);
  MenuAdd(1,"Edit",NULL);
  MenuAdd(1,"Cut",NULL);
  MenuAdd(1,"Copy",NULL);
  MenuAdd(1,"Paste",NULL);*/

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
  /*  glutDestroyWindow(TheWindow);*/
  
  MemoryDebugDump();
}
/*========================================================================*/
void MainRefreshNow(void) 
{ /* should only be called by the master thread */

  CMain *I = &Main;
  if(I->SwapFlag)
	 {
	   glutSwapBuffers();
	   I->SwapFlag=false;
	 }
  if(I->DirtyFlag)
	{
	  glutPostRedisplay();
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

  if(I->SwapFlag)
	 {
		glutSwapBuffers();
		I->SwapFlag=false;
	 }

  if(I->DirtyFlag)
	 {
	   wasDirty=true;
		glutPostRedisplay();
		I->DirtyFlag=false;
	 }
  if(!wasDirty) {
	if(false&&I->IdleFlag) { /* select to avoid racing the CPU */

	  PUnlock(cLockAPI,&_save);
	  PSleep(2000);
	  PLock(cLockAPI,&_save);
	  

	  if(I->SwapFlag)
		{
		  glutSwapBuffers();
		  I->SwapFlag=false;
		}
	  if(I->DirtyFlag) 
		{
		  glutPostRedisplay();
		  I->DirtyFlag=false;
		}
	}
  }

  PUnlock(cLockAPI,&_save);
}
/*========================================================================*/

#ifndef _PYMOL_MODULE
int main(int argc, char *argv[])
{
  myArgc=argc;
  myArgv=argv;

#else
void was_main(void) 
{
  int argc = 0;
  char *argv[1],argvv[2] = "\0";
  argv[0]=argvv;

#ifdef _DRI_WORKAROUND
  dlopen("libGL.so.1",RTLD_LAZY|RTLD_GLOBAL);
#endif

#endif  
	  
  myArgc=argc;
  myArgv=argv;
  
  glutInit(&argc, argv);
  
  glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
  
  glutInitWindowPosition(0, 175);
  glutInitWindowSize(WinX, WinY);
  
  TheWindow = glutCreateWindow("PyMol Viewer");
  
  MainInit();

  PInit();

  glutDisplayFunc(         MainDraw );
  glutReshapeFunc(         MainReshape );
  glutKeyboardFunc(        MainKey );
  glutMouseFunc(           MainButton );
  glutMotionFunc(          MainDrag );
  /*  glutPassiveMotionFunc(   MainMove );*/
  glutSpecialFunc(         MainSpecial );
  glutIdleFunc(         MainBusyIdle );

  glutPostRedisplay();

  Py_UNBLOCK_THREADS;
  
  PyMOLReady = true;

  glutMainLoop();

#ifndef _PYMOL_MODULE
  return 0;
#endif

}






