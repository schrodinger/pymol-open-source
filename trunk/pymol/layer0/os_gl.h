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

#include"os_predef.h"
#include"os_proprietary.h"

#ifndef _PYMOL_OSX

#ifdef _PYMOL_OPENGL_SHADERS
#ifndef WIN32
#define GL_GLEXT_PROTOTYPES
#endif
#include<GL/gl.h>
#include<GL/glu.h>
#include<GL/glext.h>
#else
#include<GL/gl.h>
#include<GL/glu.h>
#endif

#else

/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef _MACPYMOL_XCODE
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <OpenGL/glext.h>
#else
#include<gl.h>
#include<glu.h>
#include<glext.h>
#endif
/* END PROPRIETARY CODE SEGMENT */
#endif

/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#ifdef WIN32
#include<GL/glu.h>
#ifdef _PYMOL_OPENGL_SHADERS
#include<GL/glext.h>
#endif
#endif
/* END PROPRIETARY CODE SEGMENT */

void PyMOLReadPixels(GLint x,
                  GLint y,
                  GLsizei width,
                  GLsizei height,
                  GLenum format,
                  GLenum type,
                      GLvoid *pixels);


void PyMOLDrawPixels(GLsizei width,
                  GLsizei height,
                  GLenum format,
                  GLenum type,
                  const GLvoid *pixels);

#define P_GLUT_BUTTON_SCROLL_FORWARD  3
#define P_GLUT_BUTTON_SCROLL_BACKWARD 4
#define P_GLUT_DOUBLE_LEFT 5
#define P_GLUT_DOUBLE_MIDDLE 6
#define P_GLUT_DOUBLE_RIGHT 7
#define P_GLUT_SINGLE_LEFT 8
#define P_GLUT_SINGLE_MIDDLE 9
#define P_GLUT_SINGLE_RIGHT 10

int PyMOLCheckOpenGLErr(char *pos);

/* determine whether or not we have a real GLUT */

#ifdef _PYMOL_NO_GLUT
#define _PYMOL_PRETEND_GLUT
#endif

#ifdef _PYMOL_ACTIVEX_OLD
#define _PYMOL_WX_GLUT
#endif

/*
#ifdef _EPYMOL
#define _PYMOL_WX_GLUT
#endif
*/

#ifdef _PYMOL_MIN_GLUT
#define _PYMOL_PRETEND_GLUT
#endif

#ifdef _PYMOL_WX_GLUT
#define _PYMOL_PRETEND_GLUT
#define _PYMOL_PRETEND_GLUT_FONT
#endif

#ifdef _PYMOL_NO_MAIN
#define _PYMOL_PRETEND_GLUT
#define _PYMOL_NO_GLUT
#endif

#ifndef _PYMOL_PRETEND_GLUT

/* ============ REAL GLUT BEING USED ============= */

#ifndef _PYMOL_OSX
#include<GL/glut.h>
#else
#include<glut.h>
#endif

#ifdef FREEGLUT
#include<GL/freeglut_ext.h>
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
#define P_GLUT_MULTISAMPLE              GLUT_MULTISAMPLE
#define P_GLUT_STENCIL                  GLUT_STENCIL

#define P_GLUT_WINDOW_X                 GLUT_WINDOW_X
#define P_GLUT_WINDOW_Y                 GLUT_WINDOW_Y
#define P_GLUT_WINDOW_WIDTH             GLUT_WINDOW_WIDTH
#define P_GLUT_WINDOW_HEIGHT            GLUT_WINDOW_HEIGHT
#define P_GLUT_SCREEN_HEIGHT            GLUT_SCREEN_HEIGHT
#define P_GLUT_SCREEN_WIDTH             GLUT_SCREEN_WIDTH
#define P_GLUT_WINDOW_BORDER_WIDTH      GLUT_WINDOW_BORDER_WIDTH
#define P_GLUT_WINDOW_HEADER_HEIGHT     GLUT_WINDOW_HEADER_HEIGHT

#define p_glutGameModeString       glutGameModeString
#define p_glutEnterGameMode        glutEnterGameMode
#define p_glutLeaveGameMode        glutLeaveGameMode

#define p_glutBitmapCharacter      glutBitmapCharacter
#define p_glutSwapBuffers          glutSwapBuffers

#define p_glutCreateWindow         glutCreateWindow
#define p_glutPopWindow            glutPopWindow
#define p_glutShowWindow           glutShowWindow
#define p_glutHideWindow           glutHideWindow
#define p_glutReshapeWindow        glutReshapeWindow
#define p_glutDestroyWindow        glutDestroyWindow

#define p_glutFullScreen           glutFullScreen
#define p_glutPostRedisplay        glutPostRedisplay

#define p_glutInit                 glutInit
#define p_glutInitDisplayMode      glutInitDisplayMode
#define p_glutInitWindowPosition   glutInitWindowPosition
#define p_glutInitWindowSize       glutInitWindowSize
#define p_glutPositionWindow       glutPositionWindow

#define p_glutGet                  glutGet
#define p_glutGetModifiers         glutGetModifiers

#define p_glutDisplayFunc          glutDisplayFunc
#define p_glutReshapeFunc          glutReshapeFunc
#define p_glutKeyboardFunc         glutKeyboardFunc
#define p_glutMouseFunc            glutMouseFunc
#define p_glutMotionFunc           glutMotionFunc
#define p_glutPassiveMotionFunc    glutPassiveMotionFunc
#define p_glutSpecialFunc          glutSpecialFunc
#define p_glutIdleFunc             glutIdleFunc

#define p_glutMainLoop             glutMainLoop

#else

/* ============ GLUT EMULATION MODE ============= */

#ifndef _PYMOL_NO_GLUT

#define P_GLUT_IDLE_EVENT            0
#define P_GLUT_DISPLAY_EVENT         1
#define P_GLUT_RESHAPE_EVENT         2
#define P_GLUT_MOUSE_EVENT           3
#define P_GLUT_MOTION_EVENT          4
#define P_GLUT_CHAR_EVENT            5
#define P_GLUT_SPECIAL_EVENT         6
#define P_GLUT_PASSIVE_MOTION_EVENT  7


typedef struct {
  int event_code;
  int x,y;
  int input,state,mod;
} p_glut_event;

/* here is the pretend GLUT event handler */

void p_glutHandleEvent(p_glut_event *ev); 
int p_glutGetRedisplay(void);

#endif

/* here is the interface and constants for pretend or no GLUT */

#define P_GLUT_RGBA                     0
#define P_GLUT_DOUBLE                   2
#define P_GLUT_ACTIVE_ALT               4
#define P_GLUT_ACTIVE_CTRL              2
#define P_GLUT_ACTIVE_SHIFT             1
#define P_GLUT_BITMAP_8_BY_13           ((void*)3)
#define P_GLUT_DEPTH                    16
#define P_GLUT_DISPLAY_MODE_POSSIBLE    400

#define P_GLUT_DOWN           0
#define P_GLUT_UP             1

#define P_GLUT_STEREO         256
#define P_GLUT_MULTISAMPLE    128
#define P_GLUT_STENCIL        32

#define P_GLUT_KEY_F1         1
#define P_GLUT_KEY_F2         2
#define P_GLUT_KEY_F3         3
#define P_GLUT_KEY_F4         4
#define P_GLUT_KEY_F5         5
#define P_GLUT_KEY_F6         6
#define P_GLUT_KEY_F7         7
#define P_GLUT_KEY_F8         8
#define P_GLUT_KEY_F9         9
#define P_GLUT_KEY_F10        10
#define P_GLUT_KEY_F11        11
#define P_GLUT_KEY_F12        12
#define P_GLUT_KEY_LEFT       100
#define P_GLUT_KEY_UP         101
#define P_GLUT_KEY_RIGHT      102
#define P_GLUT_KEY_DOWN       103
#define P_GLUT_KEY_PAGE_UP    104
#define P_GLUT_KEY_PAGE_DOWN  105
#define P_GLUT_KEY_HOME       106
#define P_GLUT_KEY_END        107
#define P_GLUT_KEY_INSERT     108
#define P_GLUT_LEFT_BUTTON    0
#define P_GLUT_MIDDLE_BUTTON  1
#define P_GLUT_RIGHT_BUTTON   2

#define P_GLUT_WINDOW_X                 5
#define P_GLUT_WINDOW_Y                 6
#define P_GLUT_WINDOW_WIDTH             7
#define P_GLUT_WINDOW_HEIGHT            8
#define P_GLUT_SCREEN_WIDTH             9
#define P_GLUT_SCREEN_HEIGHT            10
#define P_GLUT_WINDOW_BORDER_WIDTH             11
#define P_GLUT_WINDOW_HEADER_HEIGHT            12


#ifndef _PYMOL_NO_GLUT
void     p_glutGameModeString(char *str);
void     p_glutEnterGameMode(void);
void     p_glutLeaveGameMode(void);

void     p_glutBitmapCharacter(void *font, int character);
void     p_glutSwapBuffers(void);

int      p_glutCreateWindow(const char *title); /* NOTE: once this function called,
                                                 * we need to have a valid GL rendering context */
void     p_glutPopWindow(void);
void     p_glutShowWindow(void);
void     p_glutReshapeWindow(int width, int height);
void     p_glutDestroyWindow(int theWindow);

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

void     p_glutPositionWindow(int x,int y);
void     p_glutMainLoop(void);

#endif

#endif

#endif
