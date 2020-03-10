/**
 * Real GLUT mode
 */

#pragma once

#ifdef _PYMOL_PRETEND_GLUT
#error "do not include directly, include os_gl.h"
#endif

#ifndef _PYMOL_OSX
#include <GL/glut.h>
#else
#import <GLUT/glut.h>
#endif

#ifdef FREEGLUT
#include <GL/freeglut_ext.h>
#endif

/* These are the only glut constants and functions that PyMOL uses ... */

#define P_GLUT_ACTIVE_ALT               GLUT_ACTIVE_ALT
#define P_GLUT_ACTIVE_CTRL              GLUT_ACTIVE_CTRL
#define P_GLUT_ACTIVE_SHIFT             GLUT_ACTIVE_SHIFT
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
