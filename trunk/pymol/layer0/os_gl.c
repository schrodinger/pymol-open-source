#include"os_gl.h"

#ifdef _PYMOL_MIN_GLUT

/* NULL GLUT: Bare Minimum GLUT implementation for getting PyMOL up in 
   preparation for porting PYMOL to a non-GLUT environment (such as ActiveX)... */

#include<GL/glut.h>

#ifndef NULL
#define NULL ((void*)0)
#endif

static void (*idleFunc)(void) = NULL;
static int WinX = 640,WinY=480;

int      p_glutCreateWindow(const char *title) {

  char **myArgv,*myArgvv[2],myArgvvv[1024] = "pymol";
  int myArgc;

  myArgc=1;
  myArgvv[0]=myArgvvv;
  myArgvv[1]=(void*)0;
  myArgv=myArgvv;
  
  glutInit(&myArgc, myArgv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE );            
  glutInitWindowPosition(0, 175);
  glutInitWindowSize(WinX, WinY);
  
  return glutCreateWindow(title);
}
void     p_glutMainLoop(void) { 
  while(1) { if(idleFunc) idleFunc(); }

}

void     p_glutBitmapCharacter(void *font, int character){}
void     p_glutSwapBuffers(void){}

void     p_glutPopWindow(void){}
void     p_glutShowWindow(void){}

void     p_glutReshapeWindow(int width, int height){ WinX=width;WinY=height;}

void     p_glutFullScreen(void) {}
void     p_glutPostRedisplay(void) {}

void     p_glutInit(int *argcp, char **argv) {}
void     p_glutInitDisplayMode(unsigned int mode) {}
void     p_glutInitDisplayString(const char *string) {}
void     p_glutInitWindowPosition(int x, int y) {}
void     p_glutInitWindowSize(int width, int height) { WinX=width;WinY=height;}

int      p_glutGet(GLenum type) {return 0;}
int      p_glutGetModifiers(void) {return 0;}

void     p_glutDisplayFunc(void (*func)(void)) {}
void     p_glutReshapeFunc(void (*func)(int width, int height)) {}
void     p_glutKeyboardFunc(void (*func)(unsigned char key, int x, int y)) {}
void     p_glutMouseFunc(void (*func)(int button, int state, int x, int y)) {}
void     p_glutMotionFunc(void (*func)(int x, int y)) {}
void     p_glutSpecialFunc(void (*func)(int key, int x, int y)) {}
void     p_glutIdleFunc(void (*func)(void)) { idleFunc = func; }

#endif


#ifdef _PYMOL_WX_GLUT

#include<GL/glut.h>
#include"P.h"
#include"main.h"

/* wxWindow-based GLUT emulation layer */

#ifndef NULL
#define NULL ((void*)0)
#endif

static void (*idleFunc)(void) = NULL;
static void (*displayFunc)(void) = NULL;
static void (*reshapeFunc)(int width,int height) = NULL;
static int WinX = 640,WinY=480;

void p_glutHandleEvent(p_glut_event *ev) {
  /* NOTE: this function should be called with an unblocked interpreter,
   just as if it were an asychronous thread */

  switch(ev->event_code) {
  case P_GLUT_IDLE_EVENT:
    if(idleFunc) idleFunc();
    break;
  case P_GLUT_DISPLAY_EVENT:
    if(displayFunc) displayFunc();
    break;
  case P_GLUT_RESHAPE_EVENT:
    if(reshapeFunc) reshapeFunc(640,480);
    break;
  }
}

int      p_glutCreateWindow(const char *title) {
  return 0;
}

void     p_glutMainLoop(void) { 
  /* TRICK: for "pretend GLUT", this routine must return, so that the
     host windowing system can keep control of the main thread...*/
}
void     p_glutBitmapCharacter(void *font, int character){}
void     p_glutSwapBuffers(void){}

void     p_glutPopWindow(void){}
void     p_glutShowWindow(void){}

void     p_glutReshapeWindow(int width, int height){ WinX=width;WinY=height;}

void     p_glutFullScreen(void) {}
void     p_glutPostRedisplay(void) {}

void     p_glutInit(int *argcp, char **argv) {}
void     p_glutInitDisplayMode(unsigned int mode) {}
void     p_glutInitDisplayString(const char *string) {}
void     p_glutInitWindowPosition(int x, int y) {}
void     p_glutInitWindowSize(int width, int height) { WinX=width;WinY=height;}

int      p_glutGet(GLenum type) {return 0;}
int      p_glutGetModifiers(void) {return 0;}

void     p_glutDisplayFunc(void (*func)(void)) { displayFunc = func;}
void     p_glutReshapeFunc(void (*func)(int width, int height)) { reshapeFunc = func;}
void     p_glutKeyboardFunc(void (*func)(unsigned char key, int x, int y)) {}
void     p_glutMouseFunc(void (*func)(int button, int state, int x, int y)) {}
void     p_glutMotionFunc(void (*func)(int x, int y)) {}
void     p_glutSpecialFunc(void (*func)(int key, int x, int y)) {}
void     p_glutIdleFunc(void (*func)(void)) { idleFunc = func; }

#endif

