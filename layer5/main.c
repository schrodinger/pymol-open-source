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

static GLint WinX = 640+cOrthoRightSceneMargin;
static GLint WinY = 480+cOrthoBottomSceneMargin;
static GLint Modifiers = 0;

static char **myArgv;
static int myArgc;

static int FinalInitFlag=1;

/*========================================================================*/
void MainTest(void)
{
}
/*========================================================================*/
static void MainButton(int button,int state,int x,int y)
{
  static int glMod;  
  y=WinY-y;

  glMod = glutGetModifiers();
  Modifiers = ((glMod&GLUT_ACTIVE_SHIFT) ? cOrthoSHIFT : 0) |
	 ((glMod&GLUT_ACTIVE_CTRL) ? cOrthoCTRL : 0) |
	 ((glMod&GLUT_ACTIVE_ALT) ? cOrthoALT : 0);

  if(!OrthoButton(button,state,x,y,Modifiers))
    {
    }

  if(ControlIdling()) {
	 glutIdleFunc(MainBusyIdle);
	 SceneRestartTimers();
  }
}
/*========================================================================*/
static void MainDrag(int x,int y)
{
  y=WinY-y;
  if(!OrthoDrag(x,y,Modifiers))
    {
	 }
}
/*========================================================================*/
static void MainMove(int x,int y)
{
  y=WinY-y;
  OrthoCursor(x,y);
}
/*========================================================================*/
static void MainDraw(void)
{
  int a,l;
  char *p;
  char buffer[300];
  OrthoBusyPrime();
  ExecutiveDrawNow();
  if(FinalInitFlag)
	 {
		FinalInitFlag=0;
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
		OrthoRestorePrompt();
	 }
  if(ControlIdling()) 
	 glutIdleFunc(MainBusyIdle);
}
/*========================================================================*/
static void MainKey(unsigned char k, int x, int y)
{
  switch (k) 
	 {
	 case 27:  /* Escape */
		PExit(EXIT_SUCCESS);
		break;
	 default:
		OrthoKey(k,x,y);
		return;
	 }
  glutPostRedisplay();
}

/*========================================================================*/
static void MainSpecial(int k, int x, int y)
{
  switch (k)
	 {
	 default:
		return;
	 }
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
  /*  glutDestroyWindow(TheWindow);*/
  
  MemoryDebugDump();
}

/*========================================================================*/
void MainBusyIdle(void) 
{

  if(!ControlIdling()) glutIdleFunc(NULL);
  SceneIdle();
}
/*========================================================================*/
int main(int argc, char *argv[])
{
  myArgc=argc;
  myArgv=argv;
  
  glutInit(&argc, argv);

  glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);

  glutInitWindowPosition(0, 0);
  glutInitWindowSize(WinX, WinY);

  TheWindow = glutCreateWindow("PyMol");

  
  MainInit();

  PInit();

  glutDisplayFunc(         MainDraw );
  glutReshapeFunc(         MainReshape );
  glutKeyboardFunc(        MainKey );
  glutMouseFunc(           MainButton );
  glutMotionFunc(          MainDrag );
  glutPassiveMotionFunc(   MainMove );
  glutSpecialFunc(         MainSpecial );
  glutIdleFunc(         MainBusyIdle );

  glutPostRedisplay();
  glutMainLoop();

  return 0;
}






