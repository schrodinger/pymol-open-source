#ifndef _H_os_gl
#define _H_os_gl

#include"os_predef.h"
#include"os_proprietary.h"

// hardcode either true, or (x)
#define ALWAYS_IMMEDIATE_OR(x) true

  #define _PYMOL_NO_AA_SHADERS

#ifndef PURE_OPENGL_ES_2
#define _PYMOL_ARB_SHADERS
#endif

#define GL_DEFAULT_SHADER_WITH_SETTINGS  0xffe0
#define GL_SPHERE_SHADER  0xffe1
#define GL_CYLINDER_SHADER  0xffe2
#define GL_TWO_SIDED_LIGHTING 0xffe3
#define GL_MESH_LIGHTING 0xffe4
#define GL_DOT_LIGHTING 0xffe5
#define GL_LABEL_FLOAT_TEXT 0xffe6
#define GL_DASH_TRANSPARENCY_DEPTH_TEST 0xffe7
#define GL_BACK_FACE_CULLING 0xffe8
#define GL_DEPTH_TEST_IF_FLOATING 0xffe9

// unused: 0xfff0, 0xfffc
#define GL_LABEL_SHADER  0xfffa
#define GL_BACKGROUND_SHADER  0xfffb
#define GL_DEFAULT_SHADER  0xfffd
#define GL_SHADER_LIGHTING 0xfffe
#define GL_SCREEN_SHADER  0xfff1
#define GL_RAMP_SHADER  0xfff2
#define GL_CONNECTOR_SHADER  0xfff3
#define GL_FXAA_SHADER  0xfff4
#define GL_SMAA1_SHADER  0xfff5
#define GL_SMAA2_SHADER  0xfff6
#define GL_SMAA3_SHADER  0xfff7
#define GL_TRILINES_SHADER 0xfff8
#define GL_OIT_SHADER  0xfff9
#define GL_OIT_COPY_SHADER  0xffea
#define GL_SURFACE_SHADER  0xffeb
#define GL_LINE_SHADER  0xffec

#define CGO_GL_LIGHTING 0xffef

#ifndef _PYMOL_OSX


#ifndef WIN32
#define GL_GLEXT_PROTOTYPES
#endif

#ifndef GLEW_NO_GLU
#define GLEW_NO_GLU
#endif

#include<GL/glew.h>
#include<GL/gl.h>

#define GLDOUBLEMULTMATRIX glMultMatrixd
#define GLDOUBLETRANSLATE glTranslated

#else


/* BEGIN PROPRIETARY CODE SEGMENT (see disclaimer in "os_proprietary.h") */
#include<GL/glew.h>
#include <OpenGL/gl.h>
#include <OpenGL/glext.h>
#define GLDOUBLEMULTMATRIX glMultMatrixd
#define GLDOUBLETRANSLATE glTranslated

/* END PROPRIETARY CODE SEGMENT */
#endif

#ifdef PURE_OPENGL_ES_2
#define GLLIGHTMODELI(arg1, arg2)  /* nothing */
#else
#define GLLIGHTMODELI glLightModeli
#endif

void PyMOLReadPixels(GLint x,
                     GLint y,
                     GLsizei width,
                     GLsizei height, GLenum format, GLenum type, GLvoid * pixels);

void PyMOLDrawPixels(GLsizei width,
                     GLsizei height, GLenum format, GLenum type, const GLvoid * pixels);

#define P_GLUT_BUTTON_SCROLL_FORWARD  3
#define P_GLUT_BUTTON_SCROLL_BACKWARD 4
#define P_GLUT_DOUBLE_LEFT 200
#define P_GLUT_DOUBLE_MIDDLE 201
#define P_GLUT_DOUBLE_RIGHT 202
#define P_GLUT_SINGLE_LEFT 100
#define P_GLUT_SINGLE_MIDDLE 101
#define P_GLUT_SINGLE_RIGHT 102

int PyMOLCheckOpenGLErr(const char *pos);


/* determine whether or not we have a real GLUT */

#ifdef _PYMOL_NO_MAIN
#define _PYMOL_NO_GLUT
#endif

#ifdef _PYMOL_NO_GLUT
#define _PYMOL_PRETEND_GLUT
#endif

#ifndef _PYMOL_PRETEND_GLUT


/* ============ REAL GLUT BEING USED ============= */

#include <stdlib.h>

#ifndef _PYMOL_OSX
#include<GL/glut.h>
#else
#import <GLUT/glut.h>
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
#define P_GLUT_ALPHA                    GLUT_ALPHA
#define P_GLUT_RIGHT_BUTTON             GLUT_RIGHT_BUTTON
#define P_GLUT_STEREO                   GLUT_STEREO
#define P_GLUT_UP                       GLUT_UP
#define P_GLUT_MULTISAMPLE              GLUT_MULTISAMPLE
#define P_GLUT_STENCIL                  GLUT_STENCIL
#define P_GLUT_ACCUM                    GLUT_ACCUM

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
/* we're pretending we have glut...*/
int p_glutGet(GLenum type);

/* ============ GLUT EMULATION MODE ============= */

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
void p_glutGameModeString(char *str);
void p_glutEnterGameMode(void);
void p_glutLeaveGameMode(void);

void p_glutBitmapCharacter(void *font, int character);
void p_glutSwapBuffers(void);

int p_glutCreateWindow(const char *title);      /* NOTE: once this function called,
                                                 * we need to have a valid GL rendering context */
void p_glutPopWindow(void);
void p_glutShowWindow(void);
void p_glutReshapeWindow(int width, int height);
void p_glutDestroyWindow(int theWindow);

void p_glutFullScreen(void);
void p_glutPostRedisplay(void);

void p_glutInit(int *argcp, char **argv);
void p_glutInitDisplayMode(unsigned int mode);
void p_glutInitDisplayString(const char *string);
void p_glutInitWindowPosition(int x, int y);
void p_glutInitWindowSize(int width, int height);

int p_glutGet(GLenum type);
int p_glutGetModifiers(void);

void p_glutDisplayFunc(void (*func) (void));
void p_glutReshapeFunc(void (*func) (int width, int height));
void p_glutKeyboardFunc(void (*func) (unsigned char key, int x, int y));
void p_glutMouseFunc(void (*func) (int button, int state, int x, int y));
void p_glutMotionFunc(void (*func) (int x, int y));
void p_glutSpecialFunc(void (*func) (int key, int x, int y));
void p_glutIdleFunc(void (*func) (void));

void p_glutPositionWindow(int x, int y);
void p_glutMainLoop(void);

#endif

#endif


#define GL_C_INT_TYPE uint
#define GL_C_INT_ENUM GL_UNSIGNED_INT
#define SceneGLClearColor(red,green,blue,alpha) glClearColor(red,green,blue,alpha);

#ifndef GL_FRAGMENT_PROGRAM_ARB
#define GL_FRAGMENT_PROGRAM_ARB                         0x8804
#endif

#ifndef GLAPIENTRY
#define GLAPIENTRY
#endif

#define hasFrameBufferBinding() false

#define GL_DEBUG_PUSH(title) \
  GLEW_KHR_debug ? glPushDebugGroup(GL_DEBUG_SOURCE_APPLICATION, 0, -1, title) : (void)0

#define GL_DEBUG_POP() \
  GLEW_KHR_debug ? glPopDebugGroup() : (void)0

#ifdef __cplusplus

class glDebugBlock {
public:
  explicit glDebugBlock(char const* title) {
    GL_DEBUG_PUSH(title);
  }
  ~glDebugBlock() {
    GL_DEBUG_POP();
  }
};

#define GL_DEBUG_FUN() \
  glDebugBlock glDebugBlockVariable(__FUNCTION__)

#endif /* __cplusplus */

#endif
