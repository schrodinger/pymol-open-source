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
#ifndef _H_os_gl
#define _H_os_gl

#ifdef WIN32
#include<windows.h>
#endif

#ifndef _PYMOL_OSX
#include<GL/gl.h>
#else
#include<gl.h>
#endif

/* determine whether or not we have a real GLUT */

#ifdef _PYMOL_NO_GLUT
#define _PYMOL_PRETEND_GLUT
#endif

#ifdef _PYMOL_ACTIVEX
#define _PYMOL_PRETEND_GLUT
#define _PYMOL_PRETEND_GLUT_FONT
#endif

#ifdef _PYMOL_MIN_GLUT
#define _PYMOL_PRETEND_GLUT
#endif

#ifdef _PYMOL_WX_GLUT
#define _PYMOL_PRETEND_GLUT
#define _PYMOL_PRETEND_GLUT_FONT
#endif

#ifndef _PYMOL_PRETEND_GLUT

/* real GLUT being used... */

#ifndef _PYMOL_OSX
#include<GL/glut.h>
#else
#include<glut.h>
#endif

/* These are the only glut constants and functions that PyMOL uses ... */

#define P_GLUT_ACTIVE_ALT               GLUT_ACTIVE_ALT                 
#define P_GLUT_ACTIVE_CTRL              GLUT_ACTIVE_CTRL                
#define P_GLUT_ACTIVE_SHIFT             GLUT_ACTIVE_SHIFT               
#define P_GLUT_BITMAP_8_BY_13           GLUT_BITMAP_8_BY_13             
#define P_GLUT_DEPTH                    GLUT_DEPTH
#define P_GLUT_DISPLAY_MODE_POSSIBLE    GLUT_DISPLAY_MODE_POSSIBLE      
#define P_GLUT_DOUBLE                   GLUT_DOUBLE                     
#define P_GLUT_DOWN                     GLUT_DOWN                       
#define P_GLUT_KEY_DOWN                 GLUT_KEY_DOWN                   
#define P_GLUT_KEY_LEFT                 GLUT_KEY_LEFT                   
#define P_GLUT_KEY_RIGHT                GLUT_KEY_RIGHT                  
#define P_GLUT_KEY_UP                   GLUT_KEY_UP                     
#define P_GLUT_LEFT_BUTTON              GLUT_LEFT_BUTTON                
#define P_GLUT_MIDDLE_BUTTON            GLUT_MIDDLE_BUTTON              
#define P_GLUT_RGBA                     GLUT_RGBA
#define P_GLUT_RIGHT_BUTTON             GLUT_RIGHT_BUTTON               
#define P_GLUT_STEREO                   GLUT_STEREO                     
#define P_GLUT_UP                       GLUT_UP                           

#define p_glutBitmapCharacter      glutBitmapCharacter
#define p_glutSwapBuffers          glutSwapBuffers

#define p_glutCreateWindow         glutCreateWindow
#define p_glutPopWindow            glutPopWindow
#define p_glutShowWindow           glutShowWindow
#define p_glutReshapeWindow        glutReshapeWindow

#define p_glutFullScreen           glutFullScreen
#define p_glutPostRedisplay        glutPostRedisplay

#define p_glutInit                 glutInit
#define p_glutInitDisplayMode      glutInitDisplayMode
#define p_glutInitWindowPosition   glutInitWindowPosition
#define p_glutInitWindowSize       glutInitWindowSize

#define p_glutGet                  glutGet
#define p_glutGetModifiers         glutGetModifiers

#define p_glutDisplayFunc          glutDisplayFunc
#define p_glutReshapeFunc          glutReshapeFunc
#define p_glutKeyboardFunc         glutKeyboardFunc
#define p_glutMouseFunc            glutMouseFunc
#define p_glutMotionFunc           glutMotionFunc
#define p_glutSpecialFunc          glutSpecialFunc
#define p_glutIdleFunc             glutIdleFunc

#define p_glutMainLoop             glutMainLoop

#else

#define P_GLUT_IDLE_EVENT            0
#define P_GLUT_DISPLAY_EVENT         1
#define P_GLUT_RESHAPE_EVENT         2
#define P_GLUT_MOUSE_EVENT           3
#define P_GLUT_MOTION_EVENT          4
#define P_GLUT_CHAR_EVENT            5

typedef struct {
  int event_code;
  int x,y;
  int input,state,mod;
} p_glut_event;

/* here is the pretend GLUT event handler */

void p_glutHandleEvent(p_glut_event *ev); 
int p_glutGetRedisplay(void);

/* here is the interface and constants for a pretend GLUT */

#define P_GLUT_ACTIVE_ALT               32
#define P_GLUT_ACTIVE_CTRL              64
#define P_GLUT_ACTIVE_SHIFT             128
#define P_GLUT_BITMAP_8_BY_13           ((void*)4)
#define P_GLUT_DEPTH                    5
#define P_GLUT_DISPLAY_MODE_POSSIBLE    6
#define P_GLUT_DOUBLE                   7
#define P_GLUT_RGBA                     8
#define P_GLUT_DOWN                     9
#define P_GLUT_UP                       10
#define P_GLUT_KEY_DOWN                 11
#define P_GLUT_KEY_LEFT                 12
#define P_GLUT_KEY_RIGHT                13
#define P_GLUT_KEY_UP                   14
#define P_GLUT_LEFT_BUTTON              15
#define P_GLUT_MIDDLE_BUTTON            16
#define P_GLUT_RIGHT_BUTTON             17
#define P_GLUT_STEREO                   18

void     p_glutBitmapCharacter(void *font, int character);
void     p_glutSwapBuffers(void);

int      p_glutCreateWindow(const char *title); /* NOTE: once this function called,
                                                 * we need to have a valid GL rendering context */
void     p_glutPopWindow(void);
void     p_glutShowWindow(void);
void     p_glutReshapeWindow(int width, int height);

void     p_glutFullScreen(void);
void     p_glutPostRedisplay(void);

void     p_glutInit(int *argcp, char **argv);
void     p_glutInitDisplayMode(unsigned int mode);
void     p_glutInitDisplayString(const char *string);
void     p_glutInitWindowPosition(int x, int y);
void     p_glutInitWindowSize(int width, int height);

int      p_glutGet(GLenum type);
int      p_glutGetModifiers(void);

void     p_glutDisplayFunc(void (*func)(void));
void     p_glutReshapeFunc(void (*func)(int width, int height));
void     p_glutKeyboardFunc(void (*func)(unsigned char key, int x, int y));
void     p_glutMouseFunc(void (*func)(int button, int state, int x, int y));
void     p_glutMotionFunc(void (*func)(int x, int y));
void     p_glutSpecialFunc(void (*func)(int key, int x, int y));
void     p_glutIdleFunc(void (*func)(void));

void     p_glutMainLoop(void);

#endif

#endif
