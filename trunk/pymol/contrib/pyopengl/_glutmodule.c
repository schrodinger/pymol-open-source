#ifdef _PYMOL_WIN32
#include"os_predef.h"
#endif

#ifndef _H_PYMOL_NO_GLUT

/**
 *
 * GLUT Module for Python
 * 
 * Version: 0.7.PyMOL (Modified for Usage with PyMOL)
 *
 * Date: 96/09/12 (September)
 *
 * Authors: Tom Schwaller     <tom.schwaller@linux-magazin.de>
 *          Jim Hugunin       <hugunin@xerox.parc.com>
 *          David Ascher      <da@ski.org>
 *          Michel Sanner     <sanner@scripps.edu>
 * 
 * Notes:
 *
 *   - the file glutconst.py is also needed, as it defines all kinds of
 *     useful constants
 *
 *   - The SGISTEREO define should be on if you want to have access
 *     to the SGI fullscreen stereo control from Python.  This
 *     defines four more calls: start_stereo, stop_stereo, left and right.
 *
 *   - This module should compile as is w/ both GLUT 2 and GLUT 3.
 *
 *   - Note that the Setup line should read something like:
 *
 *     opengl openglmodule.c glut_shapes.c -lGL -lGLU -lnumpymodule
 *     glu glumodule.c  -lGLU
 *     glut glutmodule.c -lGL -lGLU -lglut -lXi -lXmu 
 *     
***/

/* #define SGISTEREO  */


/* ========================== Module glut =========================== */

#include "Python.h"

#ifdef WIN32
#define MS_WIN32
#endif

/*
   GLUT_API_VERSION
   Define this flag to a number that determines which functions are
   to be available. For some reason there doesn't seem to be a way 
   to check for the GLUT version actually available from the library,
   and 3.7.2 defaults to 3.6's version ID (3).

   To get an older version of the API, use GLUT_API_VERSION appropriate 
   to that version.

   See glut.h for flag descriptions.
 */

#ifndef _PYMOL_OSX
#include "GL/glut.h"
#else
#include <glut.h>
#endif


/* #include "GL/glutint.h" */

static PyObject *glut_Error;
void DeviceButtonPressGrab(void *x, void *y, void *z);
void DeviceButtonPressGrab(void *x, void *y, void *z)
{
    printf("help, shouldn't be here\n");
}

static void *map_lookup(char *mp[][2], char *nm)
{
    int i = 0;
    while (mp[i][0] != NULL) {
	if (strcmp(mp[i][0], nm) == 0) {
	    return (void *) mp[i][1];
	}
	i = i + 1;
    }
    return NULL;
}

char *glutFonts[][2] =
{
    {"glut9by15", (char *) GLUT_BITMAP_9_BY_15},
    {"glut8by13", (char *) GLUT_BITMAP_8_BY_13},
    {"glutTimesRoman10", (char *) GLUT_BITMAP_TIMES_ROMAN_10},
    {"glutTimesRoman24", (char *) GLUT_BITMAP_TIMES_ROMAN_24},
#if (GLUT_API_VERSION >= 3)
    {"glutHelvetica10", (char *) GLUT_BITMAP_HELVETICA_10},
    {"glutHelvetica12", (char *) GLUT_BITMAP_HELVETICA_12},
    {"glutHelvetica18", (char *) GLUT_BITMAP_HELVETICA_18},
    {"glutStrokeRoman", (char *) GLUT_STROKE_ROMAN},
    {"glutStrokeRomanFixed", (char *) GLUT_STROKE_MONO_ROMAN},
#endif
    {NULL, NULL},
};

static PyObject *glut_glutCreateMenuObject;

static PyObject *glut_SetglutCreateMenuCallback(PyObject * self, PyObject * args)
{
    PyObject *func;
    if (!PyArg_ParseTuple(args, "O", &func))
	return NULL;
    Py_XDECREF(glut_glutCreateMenuObject);
    glut_glutCreateMenuObject = func;
    Py_INCREF(func);
    Py_INCREF(Py_None);
    return Py_None;
}


static void glutCreateMenuCallback(int x)
{
    PyObject *_args = NULL;
    PyObject *_res = NULL;
    _args = Py_BuildValue("(i)",
			  x);
    _res = PyObject_CallObject(glut_glutCreateMenuObject, _args);
    if (_res == NULL) {
	PyErr_Print();
    }
    Py_XDECREF(_res);
}

static PyObject *glut_glutDisplayFuncObject;

static PyObject *glut_SetglutDisplayFuncCallback(PyObject * self, PyObject * args)
{
    PyObject *func;
    if (!PyArg_ParseTuple(args, "O", &func))
	return NULL;
    Py_XDECREF(glut_glutDisplayFuncObject);
    glut_glutDisplayFuncObject = func;
    Py_INCREF(func);
    Py_INCREF(Py_None);
    return Py_None;
}


static void glutDisplayFuncCallback(void)
{
    PyObject *_args = NULL;
    PyObject *_res = NULL;
    _args = Py_BuildValue("()");
    _res = PyObject_CallObject(glut_glutDisplayFuncObject, _args);
    if (_res == NULL) {
	PyErr_Print();
    }
    Py_XDECREF(_res);
}

static PyObject *glut_glutReshapeFuncObject;

static PyObject *glut_SetglutReshapeFuncCallback(PyObject * self, PyObject * args)
{
    PyObject *func;
    if (!PyArg_ParseTuple(args, "O", &func))
	return NULL;
    Py_XDECREF(glut_glutReshapeFuncObject);
    glut_glutReshapeFuncObject = func;
    Py_INCREF(func);
    Py_INCREF(Py_None);
    return Py_None;
}


static void glutReshapeFuncCallback(int width,
				    int height)
{
    PyObject *_args = NULL;
    PyObject *_res = NULL;
    _args = Py_BuildValue("(ii)",
			  width,
			  height);
    _res = PyObject_CallObject(glut_glutReshapeFuncObject, _args);
    if (_res == NULL) {
	PyErr_Print();
    }
    Py_XDECREF(_res);
}

static PyObject *glut_glutKeyboardFuncObject;

static PyObject *glut_SetglutKeyboardFuncCallback(PyObject * self, PyObject * args)
{
    PyObject *func;
    if (!PyArg_ParseTuple(args, "O", &func))
	return NULL;
    Py_XDECREF(glut_glutKeyboardFuncObject);
    glut_glutKeyboardFuncObject = func;
    Py_INCREF(func);
    Py_INCREF(Py_None);
    return Py_None;
}

#if (GLUT_API_VERSION >= 4)
static PyObject *glut_glutKeyboardUpFuncObject;

static PyObject *glut_SetglutKeyboardUpFuncCallback(PyObject * self, PyObject * args)
{
    PyObject *func;
    if (!PyArg_ParseTuple(args, "O", &func))
	return NULL;
    Py_XDECREF(glut_glutKeyboardUpFuncObject);
    glut_glutKeyboardUpFuncObject = func;
    Py_INCREF(func);
    Py_INCREF(Py_None);
    return Py_None;
}
#endif

static void glutKeyboardFuncCallback(unsigned char key,
				     int x,
				     int y)
{
    PyObject *_args = NULL;
    PyObject *_res = NULL;
    _args = Py_BuildValue("(cii)",
			  key,
			  x,
			  y);
    _res = PyObject_CallObject(glut_glutKeyboardFuncObject, _args);
    if (_res == NULL) {
	PyErr_Print();
    }
    Py_XDECREF(_res);
}

#if (GLUT_API_VERSION >= 4)
static void glutKeyboardUpFuncCallback(unsigned char key,
				       int x,
				       int y)
{
    PyObject *_args = NULL;
    PyObject *_res = NULL;
    _args = Py_BuildValue("(cii)",
			  key,
			  x,
			  y);
    _res = PyObject_CallObject(glut_glutKeyboardUpFuncObject, _args);
    if (_res == NULL) {
	PyErr_Print();
    }
    Py_XDECREF(_res);
}

static PyObject *glut_glutIgnoreKeyRepeat(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    int newvalue;
    if (!PyArg_ParseTuple(_args, "i", &newvalue))
	return NULL;
    glutIgnoreKeyRepeat(newvalue);
    printf("ignore key repeat:%d\n", newvalue);
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}
#endif

static PyObject *glut_glutMouseFuncObject;

static PyObject *glut_SetglutMouseFuncCallback(PyObject * self, PyObject * args)
{
    PyObject *func;
    if (!PyArg_ParseTuple(args, "O", &func))
	return NULL;
    Py_XDECREF(glut_glutMouseFuncObject);
    glut_glutMouseFuncObject = func;
    Py_INCREF(func);
    Py_INCREF(Py_None);
    return Py_None;
}


static void glutMouseFuncCallback(int button,
				  int state,
				  int x,
				  int y)
{
    PyObject *_args = NULL;
    PyObject *_res = NULL;
    _args = Py_BuildValue("(iiii)",
			  button,
			  state,
			  x,
			  y);
    _res = PyObject_CallObject(glut_glutMouseFuncObject, _args);
    if (_res == NULL) {
	PyErr_Print();
    }
    Py_XDECREF(_res);
}

static PyObject *glut_glutMotionFuncObject;

static PyObject *glut_SetglutMotionFuncCallback(PyObject * self, PyObject * args)
{
    PyObject *func;
    if (!PyArg_ParseTuple(args, "O", &func))
	return NULL;
    Py_XDECREF(glut_glutMotionFuncObject);
    glut_glutMotionFuncObject = func;
    Py_INCREF(func);
    Py_INCREF(Py_None);
    return Py_None;
}


static void glutMotionFuncCallback(int x,
				   int y)
{
    PyObject *_args = NULL;
    PyObject *_res = NULL;
    _args = Py_BuildValue("(ii)",
			  x,
			  y);
    _res = PyObject_CallObject(glut_glutMotionFuncObject, _args);
    if (_res == NULL) {
	PyErr_Print();
    }
    Py_XDECREF(_res);
}

static PyObject *glut_glutPassiveMotionFuncObject;

static PyObject *glut_SetglutPassiveMotionFuncCallback(PyObject * self, PyObject * args)
{
    PyObject *func;
    if (!PyArg_ParseTuple(args, "O", &func))
	return NULL;
    Py_XDECREF(glut_glutPassiveMotionFuncObject);
    glut_glutPassiveMotionFuncObject = func;
    Py_INCREF(func);
    Py_INCREF(Py_None);
    return Py_None;
}


static void glutPassiveMotionFuncCallback(int x,
					  int y)
{
    PyObject *_args = NULL;
    PyObject *_res = NULL;
    _args = Py_BuildValue("(ii)",
			  x,
			  y);
    _res = PyObject_CallObject(glut_glutPassiveMotionFuncObject, _args);
    if (_res == NULL) {
	PyErr_Print();
    }
    Py_XDECREF(_res);
}

static PyObject *glut_glutEntryFuncObject;

static PyObject *glut_SetglutEntryFuncCallback(PyObject * self, PyObject * args)
{
    PyObject *func;
    if (!PyArg_ParseTuple(args, "O", &func))
	return NULL;
    Py_XDECREF(glut_glutEntryFuncObject);
    glut_glutEntryFuncObject = func;
    Py_INCREF(func);
    Py_INCREF(Py_None);
    return Py_None;
}


static void glutEntryFuncCallback(int state)
{
    PyObject *_args = NULL;
    PyObject *_res = NULL;
    _args = Py_BuildValue("(i)",
			  state);
    _res = PyObject_CallObject(glut_glutEntryFuncObject, _args);
    if (_res == NULL) {
	PyErr_Print();
    }
    Py_XDECREF(_res);
}

static PyObject *glut_glutVisibilityFuncObject;

static PyObject *glut_SetglutVisibilityFuncCallback(PyObject * self, PyObject * args)
{
    PyObject *func;
    if (!PyArg_ParseTuple(args, "O", &func))
	return NULL;
    Py_XDECREF(glut_glutVisibilityFuncObject);
    glut_glutVisibilityFuncObject = func;
    Py_INCREF(func);
    Py_INCREF(Py_None);
    return Py_None;
}


static void glutVisibilityFuncCallback(int state)
{
    PyObject *_args = NULL;
    PyObject *_res = NULL;
    _args = Py_BuildValue("(i)",
			  state);
    _res = PyObject_CallObject(glut_glutVisibilityFuncObject, _args);
    if (_res == NULL) {
	PyErr_Print();
    }
    Py_XDECREF(_res);
}

static PyObject *glut_glutIdleFuncObject;

static PyObject *glut_SetglutIdleFuncCallback(PyObject * self, PyObject * args)
{
    PyObject *func;
    if (!PyArg_ParseTuple(args, "O", &func))
	return NULL;
    Py_XDECREF(glut_glutIdleFuncObject);
    glut_glutIdleFuncObject = func;
    Py_INCREF(func);
    Py_INCREF(Py_None);
    return Py_None;
}


static void glutIdleFuncCallback(void)
{
    PyObject *_args = NULL;
    PyObject *_res = NULL;
    _args = Py_BuildValue("()");
    _res = PyObject_CallObject(glut_glutIdleFuncObject, _args);
    if (_res == NULL) {
	PyErr_Print();
    }
    Py_XDECREF(_res);
}

static PyObject *glut_glutTimerFuncObject;

static PyObject *glut_SetglutTimerFuncCallback(PyObject * self, PyObject * args)
{
    PyObject *func;
    if (!PyArg_ParseTuple(args, "O", &func))
	return NULL;
    Py_XDECREF(glut_glutTimerFuncObject);
    glut_glutTimerFuncObject = func;
    Py_INCREF(func);
    Py_INCREF(Py_None);
    return Py_None;
}


static void glutTimerFuncCallback(int value)
{
    PyObject *_args = NULL;
    PyObject *_res = NULL;
    _args = Py_BuildValue("(i)",
			  value);
    _res = PyObject_CallObject(glut_glutTimerFuncObject, _args);
    if (_res == NULL) {
	PyErr_Print();
    }
    Py_XDECREF(_res);
}

static PyObject *glut_glutMenuStateFuncObject;

static PyObject *glut_SetglutMenuStateFuncCallback(PyObject * self, PyObject * args)
{
    PyObject *func;
    if (!PyArg_ParseTuple(args, "O", &func))
	return NULL;
    Py_XDECREF(glut_glutMenuStateFuncObject);
    glut_glutMenuStateFuncObject = func;
    Py_INCREF(func);
    Py_INCREF(Py_None);
    return Py_None;
}


static void glutMenuStateFuncCallback(int state)
{
    PyObject *_args = NULL;
    PyObject *_res = NULL;
    _args = Py_BuildValue("(i)",
			  state);
    _res = PyObject_CallObject(glut_glutMenuStateFuncObject, _args);
    if (_res == NULL) {
	PyErr_Print();
    }
    Py_XDECREF(_res);
}

static PyObject *glut_glutSpecialFuncObject;

static PyObject *glut_SetglutSpecialFuncCallback(PyObject * self, PyObject * args)
{
    PyObject *func;
    if (!PyArg_ParseTuple(args, "O", &func))
	return NULL;
    Py_XDECREF(glut_glutSpecialFuncObject);
    glut_glutSpecialFuncObject = func;
    Py_INCREF(func);
    Py_INCREF(Py_None);
    return Py_None;
}


static void glutSpecialFuncCallback(int key,
				    int x,
				    int y)
{
    PyObject *_args = NULL;
    PyObject *_res = NULL;
    _args = Py_BuildValue("(iii)",
			  key,
			  x,
			  y);
    _res = PyObject_CallObject(glut_glutSpecialFuncObject, _args);
    if (_res == NULL) {
	PyErr_Print();
    }
    Py_XDECREF(_res);
}

#if (GLUT_API_VERSION >= 4)
static PyObject *glut_glutSpecialUpFuncObject;

static PyObject *glut_SetglutSpecialUpFuncCallback(PyObject * self, PyObject * args)
{
    PyObject *func;
    if (!PyArg_ParseTuple(args, "O", &func))
	return NULL;
    Py_XDECREF(glut_glutSpecialUpFuncObject);
    glut_glutSpecialUpFuncObject = func;
    Py_INCREF(func);
    Py_INCREF(Py_None);
    return Py_None;
}


static void glutSpecialUpFuncCallback(int key,
				      int x,
				      int y)
{
    PyObject *_args = NULL;
    PyObject *_res = NULL;
    _args = Py_BuildValue("(iii)",
			  key,
			  x,
			  y);
    _res = PyObject_CallObject(glut_glutSpecialUpFuncObject, _args);
    if (_res == NULL) {
	PyErr_Print();
    }
    Py_XDECREF(_res);
}
#endif

static PyObject *glut_glutSpaceballMotionFuncObject;

static PyObject *glut_SetglutSpaceballMotionFuncCallback(PyObject * self, PyObject * args)
{
    PyObject *func;
    if (!PyArg_ParseTuple(args, "O", &func))
	return NULL;
    Py_XDECREF(glut_glutSpaceballMotionFuncObject);
    glut_glutSpaceballMotionFuncObject = func;
    Py_INCREF(func);
    Py_INCREF(Py_None);
    return Py_None;
}


static void glutSpaceballMotionFuncCallback(int x,
					    int y,
					    int z)
{
    PyObject *_args = NULL;
    PyObject *_res = NULL;
    _args = Py_BuildValue("(iii)",
			  x,
			  y,
			  z);
    _res = PyObject_CallObject(glut_glutSpaceballMotionFuncObject, _args);
    if (_res == NULL) {
	PyErr_Print();
    }
    Py_XDECREF(_res);
}

static PyObject *glut_glutSpaceballRotateFuncObject;

static PyObject *glut_SetglutSpaceballRotateFuncCallback(PyObject * self, PyObject * args)
{
    PyObject *func;
    if (!PyArg_ParseTuple(args, "O", &func))
	return NULL;
    Py_XDECREF(glut_glutSpaceballRotateFuncObject);
    glut_glutSpaceballRotateFuncObject = func;
    Py_INCREF(func);
    Py_INCREF(Py_None);
    return Py_None;
}


static void glutSpaceballRotateFuncCallback(int x,
					    int y,
					    int z)
{
    PyObject *_args = NULL;
    PyObject *_res = NULL;
    _args = Py_BuildValue("(iii)",
			  x,
			  y,
			  z);
    _res = PyObject_CallObject(glut_glutSpaceballRotateFuncObject, _args);
    if (_res == NULL) {
	PyErr_Print();
    }
    Py_XDECREF(_res);
}

static PyObject *glut_glutSpaceballButtonFuncObject;

static PyObject *glut_SetglutSpaceballButtonFuncCallback(PyObject * self, PyObject * args)
{
    PyObject *func;
    if (!PyArg_ParseTuple(args, "O", &func))
	return NULL;
    Py_XDECREF(glut_glutSpaceballButtonFuncObject);
    glut_glutSpaceballButtonFuncObject = func;
    Py_INCREF(func);
    Py_INCREF(Py_None);
    return Py_None;
}


static void glutSpaceballButtonFuncCallback(int button,
					    int state)
{
    PyObject *_args = NULL;
    PyObject *_res = NULL;
    _args = Py_BuildValue("(ii)",
			  button,
			  state);
    _res = PyObject_CallObject(glut_glutSpaceballButtonFuncObject, _args);
    if (_res == NULL) {
	PyErr_Print();
    }
    Py_XDECREF(_res);
}

static PyObject *glut_glutButtonBoxFuncObject;

static PyObject *glut_SetglutButtonBoxFuncCallback(PyObject * self, PyObject * args)
{
    PyObject *func;
    if (!PyArg_ParseTuple(args, "O", &func))
	return NULL;
    Py_XDECREF(glut_glutButtonBoxFuncObject);
    glut_glutButtonBoxFuncObject = func;
    Py_INCREF(func);
    Py_INCREF(Py_None);
    return Py_None;
}


static void glutButtonBoxFuncCallback(int button,
				      int state)
{
    PyObject *_args = NULL;
    PyObject *_res = NULL;
    _args = Py_BuildValue("(ii)",
			  button,
			  state);
    _res = PyObject_CallObject(glut_glutButtonBoxFuncObject, _args);
    if (_res == NULL) {
	PyErr_Print();
    }
    Py_XDECREF(_res);
}

static PyObject *glut_glutDialsFuncObject;

static PyObject *glut_SetglutDialsFuncCallback(PyObject * self, PyObject * args)
{
    PyObject *func;
    if (!PyArg_ParseTuple(args, "O", &func))
	return NULL;
    Py_XDECREF(glut_glutDialsFuncObject);
    glut_glutDialsFuncObject = func;
    Py_INCREF(func);
    Py_INCREF(Py_None);
    return Py_None;
}


static void glutDialsFuncCallback(int dial,
				  int value)
{
    PyObject *_args = NULL;
    PyObject *_res = NULL;
    _args = Py_BuildValue("(ii)",
			  dial,
			  value);
    _res = PyObject_CallObject(glut_glutDialsFuncObject, _args);
    if (_res == NULL) {
	PyErr_Print();
    }
    Py_XDECREF(_res);
}

static PyObject *glut_glutTabletMotionFuncObject;

static PyObject *glut_SetglutTabletMotionFuncCallback(PyObject * self, PyObject * args)
{
    PyObject *func;
    if (!PyArg_ParseTuple(args, "O", &func))
	return NULL;
    Py_XDECREF(glut_glutTabletMotionFuncObject);
    glut_glutTabletMotionFuncObject = func;
    Py_INCREF(func);
    Py_INCREF(Py_None);
    return Py_None;
}


static void glutTabletMotionFuncCallback(int x,
					 int y)
{
    PyObject *_args = NULL;
    PyObject *_res = NULL;
    _args = Py_BuildValue("(ii)",
			  x,
			  y);
    _res = PyObject_CallObject(glut_glutTabletMotionFuncObject, _args);
    if (_res == NULL) {
	PyErr_Print();
    }
    Py_XDECREF(_res);
}

static PyObject *glut_glutTabletButtonFuncObject;

static PyObject *glut_SetglutTabletButtonFuncCallback(PyObject * self, PyObject * args)
{
    PyObject *func;
    if (!PyArg_ParseTuple(args, "O", &func))
	return NULL;
    Py_XDECREF(glut_glutTabletButtonFuncObject);
    glut_glutTabletButtonFuncObject = func;
    Py_INCREF(func);
    Py_INCREF(Py_None);
    return Py_None;
}


static void glutTabletButtonFuncCallback(int button,
					 int state,
					 int x,
					 int y)
{
    PyObject *_args = NULL;
    PyObject *_res = NULL;
    _args = Py_BuildValue("(iiii)",
			  button,
			  state,
			  x,
			  y);
    _res = PyObject_CallObject(glut_glutTabletButtonFuncObject, _args);
    if (_res == NULL) {
	PyErr_Print();
    }
    Py_XDECREF(_res);
}

#if (GLUT_API_VERSION >= 3)
static PyObject *glut_glutMenuStatusFuncObject;

static PyObject *glut_SetglutMenuStatusFuncCallback(PyObject * self, PyObject * args)
{
    PyObject *func;
    if (!PyArg_ParseTuple(args, "O", &func))
	return NULL;
    Py_XDECREF(glut_glutMenuStatusFuncObject);
    glut_glutMenuStatusFuncObject = func;
    Py_INCREF(func);
    Py_INCREF(Py_None);
    return Py_None;
}


static void glutMenuStatusFuncCallback(int status,
				       int x,
				       int y)
{
    PyObject *_args = NULL;
    PyObject *_res = NULL;
    _args = Py_BuildValue("(iii)",
			  status,
			  x,
			  y);
    _res = PyObject_CallObject(glut_glutMenuStatusFuncObject, _args);
    if (_res == NULL) {
	PyErr_Print();
    }
    Py_XDECREF(_res);
}

static PyObject *glut_glutOverlayDisplayFuncObject;

static PyObject *glut_SetglutOverlayDisplayFuncCallback(PyObject * self, PyObject * args)
{
    PyObject *func;
    if (!PyArg_ParseTuple(args, "O", &func))
	return NULL;
    Py_XDECREF(glut_glutOverlayDisplayFuncObject);
    glut_glutOverlayDisplayFuncObject = func;
    Py_INCREF(func);
    Py_INCREF(Py_None);
    return Py_None;
}


static void glutOverlayDisplayFuncCallback(void)
{
    PyObject *_args = NULL;
    PyObject *_res = NULL;
    _args = Py_BuildValue("()");
    _res = PyObject_CallObject(glut_glutOverlayDisplayFuncObject, _args);
    if (_res == NULL) {
	PyErr_Print();
    }
    Py_XDECREF(_res);
}
#endif

static PyObject *glut_glutInit(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    int argcp = 1;
    char *argv;
#ifdef macintosh
    char **argvp = &argv;
#endif        
    if (!PyArg_ParseTuple(_args, "s", &argv))
	return NULL;
#ifdef macintosh        
    glutInit(&argcp, argvp);
#else
    glutInit(&argcp, &argv);
#endif    
    _res = Py_BuildValue("i", argcp);
    return _res;
}

static PyObject *glut_glutInitDisplayMode(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    int mode;
    if (!PyArg_ParseTuple(_args, "i",
			  &mode))
	return NULL;
    glutInitDisplayMode(mode);
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutInitWindowPosition(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    int x;
    int y;
    if (!PyArg_ParseTuple(_args, "ii",
			  &x,
			  &y))
	return NULL;
    glutInitWindowPosition(x,
			   y);
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutInitWindowSize(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    int width;
    int height;
    if (!PyArg_ParseTuple(_args, "ii",
			  &width,
			  &height))
	return NULL;
    glutInitWindowSize(width,
		       height);
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutMainLoop(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutMainLoop();
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutCreateWindow(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    int _rv;
    char *title;
    if (!PyArg_ParseTuple(_args, "s",
			  &title))
	return NULL;
    _rv = glutCreateWindow(title);
    _res = Py_BuildValue("i",
			 _rv);
    return _res;
}

static PyObject *glut_glutCreateSubWindow(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    int _rv;
    int win;
    int x;
    int y;
    int width;
    int height;
    if (!PyArg_ParseTuple(_args, "iiiii",
			  &win,
			  &x,
			  &y,
			  &width,
			  &height))
	return NULL;
    _rv = glutCreateSubWindow(win,
			      x,
			      y,
			      width,
			      height);
    _res = Py_BuildValue("i",
			 _rv);
    return _res;
}

static PyObject *glut_glutDestroyWindow(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    int win;
    if (!PyArg_ParseTuple(_args, "i",
			  &win))
	return NULL;
    glutDestroyWindow(win);
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutPostRedisplay(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutPostRedisplay();
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutSwapBuffers(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutSwapBuffers();
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutGetWindow(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    int _rv;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    _rv = glutGetWindow();
    _res = Py_BuildValue("i",
			 _rv);
    return _res;
}

static PyObject *glut_glutSetWindow(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    int win;
    if (!PyArg_ParseTuple(_args, "i",
			  &win))
	return NULL;
    glutSetWindow(win);
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutSetWindowTitle(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    char *title;
    if (!PyArg_ParseTuple(_args, "s",
			  &title))
	return NULL;
    glutSetWindowTitle(title);
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutSetIconTitle(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    char *title;
    if (!PyArg_ParseTuple(_args, "s",
			  &title))
	return NULL;
    glutSetIconTitle(title);
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutPositionWindow(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    int x;
    int y;
    if (!PyArg_ParseTuple(_args, "ii",
			  &x,
			  &y))
	return NULL;
    glutPositionWindow(x,
		       y);
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutReshapeWindow(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    int width;
    int height;
    if (!PyArg_ParseTuple(_args, "ii",
			  &width,
			  &height))
	return NULL;
    glutReshapeWindow(width,
		      height);
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutPopWindow(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutPopWindow();
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutPushWindow(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutPushWindow();
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutIconifyWindow(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutIconifyWindow();
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutShowWindow(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutShowWindow();
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutHideWindow(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutHideWindow();
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutFullScreen(PyObject * _self, PyObject * _args)
{
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutFullScreen();
    Py_INCREF(Py_None);
    return Py_None;
}
/* game mode, added in version 4 of the API */
#if (GLUT_API_VERSION >= 4)
static PyObject *glut_glutGameModeString(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    char *s;
    if (!PyArg_ParseTuple(_args, "s", &s))
	return NULL;
    glutGameModeString(s);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *glut_glutEnterGameMode(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    fprintf(stderr, "glutEnterGameMode()\n");
    glutEnterGameMode();
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *glut_glutLeaveGameMode(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    fprintf(stderr, "glutLeaveGameMode()\n");
    glutLeaveGameMode();
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *glut_glutGameModeGet(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    int _rv;
    long type;
    if (!PyArg_ParseTuple(_args, "l",
			  &type))
	return NULL;
    _rv = glutGameModeGet(type);
    _res = Py_BuildValue("i",
			 _rv);
    return _res;
}
#endif

static PyObject *glut_glutSetCursor(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    int cursor;
    if (!PyArg_ParseTuple(_args, "i",
			  &cursor))
	return NULL;
    glutSetCursor(cursor);
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutEstablishOverlay(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutEstablishOverlay();
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutRemoveOverlay(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutRemoveOverlay();
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutUseLayer(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    long layer;
    if (!PyArg_ParseTuple(_args, "l",
			  &layer))
	return NULL;
    glutUseLayer(layer);
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutPostOverlayRedisplay(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutPostOverlayRedisplay();
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutShowOverlay(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutShowOverlay();
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutHideOverlay(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutHideOverlay();
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutCreateMenu(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    int _rv;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    _rv = glutCreateMenu(glutCreateMenuCallback);
    _res = Py_BuildValue("i",
			 _rv);
    return _res;
}

static PyObject *glut_glutDestroyMenu(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    int menu;
    if (!PyArg_ParseTuple(_args, "i",
			  &menu))
	return NULL;
    glutDestroyMenu(menu);
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutGetMenu(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    int _rv;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    _rv = glutGetMenu();
    _res = Py_BuildValue("i",
			 _rv);
    return _res;
}

static PyObject *glut_glutSetMenu(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    int menu;
    if (!PyArg_ParseTuple(_args, "i",
			  &menu))
	return NULL;
    glutSetMenu(menu);
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutAddMenuEntry(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    char *label;
    int value;
    if (!PyArg_ParseTuple(_args, "si",
			  &label,
			  &value))
	return NULL;
    glutAddMenuEntry(label,
		     value);
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutAddSubMenu(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    char *label;
    int submenu;
    if (!PyArg_ParseTuple(_args, "si",
			  &label,
			  &submenu))
	return NULL;
    glutAddSubMenu(label,
		   submenu);
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutChangeToMenuEntry(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    int item;
    char *label;
    int value;
    if (!PyArg_ParseTuple(_args, "isi",
			  &item,
			  &label,
			  &value))
	return NULL;
    glutChangeToMenuEntry(item,
			  label,
			  value);
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutChangeToSubMenu(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    int item;
    char *label;
    int submenu;
    if (!PyArg_ParseTuple(_args, "isi",
			  &item,
			  &label,
			  &submenu))
	return NULL;
    glutChangeToSubMenu(item,
			label,
			submenu);
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutRemoveMenuItem(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    int item;
    if (!PyArg_ParseTuple(_args, "i",
			  &item))
	return NULL;
    glutRemoveMenuItem(item);
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutAttachMenu(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    int button;
    if (!PyArg_ParseTuple(_args, "i",
			  &button))
	return NULL;
    glutAttachMenu(button);
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutDetachMenu(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    int button;
    if (!PyArg_ParseTuple(_args, "i",
			  &button))
	return NULL;
    glutDetachMenu(button);
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

#define NULL_IF_NONE(_py_callback, _c_callback) \
                               ( _py_callback == Py_None ? NULL : _c_callback)

static PyObject *glut_glutDisplayFunc(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutDisplayFunc(NULL_IF_NONE(glut_glutDisplayFuncObject,
				 glutDisplayFuncCallback));
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutReshapeFunc(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutReshapeFunc(NULL_IF_NONE(glut_glutReshapeFuncObject,
				 glutReshapeFuncCallback));
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutKeyboardFunc(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutKeyboardFunc(NULL_IF_NONE(glut_glutKeyboardFuncObject,
				  glutKeyboardFuncCallback));
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

/* added in version 4 of the API */
#if (GLUT_API_VERSION >= 4)
static PyObject *glut_glutKeyboardUpFunc(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutKeyboardFunc(NULL_IF_NONE(glut_glutKeyboardUpFuncObject,
				  glutKeyboardUpFuncCallback));
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}
#endif

static PyObject *glut_glutMouseFunc(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutMouseFunc(NULL_IF_NONE(glut_glutMouseFuncObject,
			       glutMouseFuncCallback));
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutMotionFunc(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutMotionFunc(NULL_IF_NONE(glut_glutMotionFuncObject,
				glutMotionFuncCallback));
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutPassiveMotionFunc(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutPassiveMotionFunc(NULL_IF_NONE(glut_glutPassiveMotionFuncObject,
				       glutPassiveMotionFuncCallback));
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutEntryFunc(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutEntryFunc(NULL_IF_NONE(glut_glutEntryFuncObject,
			       glutEntryFuncCallback));
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutVisibilityFunc(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutVisibilityFunc(NULL_IF_NONE(glut_glutVisibilityFuncObject,
				    glutVisibilityFuncCallback));
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutIdleFunc(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutIdleFunc(NULL_IF_NONE(glut_glutIdleFuncObject,
			      glutIdleFuncCallback));
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutTimerFunc(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    int millis;
    int value;
    if (!PyArg_ParseTuple(_args, "ii",
			  &millis,
			  &value))
	return NULL;
    glutTimerFunc(millis,
		  NULL_IF_NONE(glut_glutTimerFuncObject,
			       glutTimerFuncCallback),
		  value);
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutMenuStateFunc(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutMenuStateFunc(NULL_IF_NONE(glut_glutMenuStateFuncObject,
				   glutMenuStateFuncCallback));
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutSpecialFunc(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutSpecialFunc(NULL_IF_NONE(glut_glutSpecialFuncObject,
				 glutSpecialFuncCallback));
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

/* added in version 4 of the API */
#if (GLUT_API_VERSION >= 4)
static PyObject *glut_glutSpecialUpFunc(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutSpecialUpFunc(NULL_IF_NONE(glut_glutSpecialUpFuncObject,
				   glutSpecialUpFuncCallback));
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}
#endif

static PyObject *glut_glutSpaceballMotionFunc(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutSpaceballMotionFunc(NULL_IF_NONE(glut_glutSpaceballMotionFuncObject,
				       glutSpaceballMotionFuncCallback));
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutSpaceballRotateFunc(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutSpaceballRotateFunc(NULL_IF_NONE(glut_glutSpaceballRotateFuncObject,
				       glutSpaceballRotateFuncCallback));
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutSpaceballButtonFunc(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutSpaceballButtonFunc(NULL_IF_NONE(glut_glutSpaceballButtonFuncObject,
				       glutSpaceballButtonFuncCallback));
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutButtonBoxFunc(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutButtonBoxFunc(NULL_IF_NONE(glut_glutButtonBoxFuncObject,
				   glutButtonBoxFuncCallback));
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutDialsFunc(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutDialsFunc(NULL_IF_NONE(glut_glutDialsFuncObject,
			       glutDialsFuncCallback));
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutTabletMotionFunc(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutTabletMotionFunc(NULL_IF_NONE(glut_glutTabletMotionFuncObject,
				      glutTabletMotionFuncCallback));
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutTabletButtonFunc(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutTabletButtonFunc(NULL_IF_NONE(glut_glutTabletButtonFuncObject,
				      glutTabletButtonFuncCallback));
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutMenuStatusFunc(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutMenuStatusFunc(NULL_IF_NONE(glut_glutMenuStatusFuncObject,
				    glutMenuStatusFuncCallback));
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutOverlayDisplayFunc(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutOverlayDisplayFunc(NULL_IF_NONE(glut_glutOverlayDisplayFuncObject,
					glutOverlayDisplayFuncCallback));
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutSetColor(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    int x;
    float red;
    float green;
    float blue;
    if (!PyArg_ParseTuple(_args, "ifff",
			  &x,
			  &red,
			  &green,
			  &blue))
	return NULL;
    glutSetColor(x,
		 red,
		 green,
		 blue);
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutGetColor(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    float _rv;
    int ndx;
    int component;
    if (!PyArg_ParseTuple(_args, "ii",
			  &ndx,
			  &component))
	return NULL;
    _rv = glutGetColor(ndx,
		       component);
    _res = Py_BuildValue("f",
			 _rv);
    return _res;
}

static PyObject *glut_glutCopyColormap(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    int win;
    if (!PyArg_ParseTuple(_args, "i",
			  &win))
	return NULL;
    glutCopyColormap(win);
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutGet(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    int _rv;
    long type;
    if (!PyArg_ParseTuple(_args, "l",
			  &type))
	return NULL;
    _rv = glutGet(type);
    _res = Py_BuildValue("i",
			 _rv);
    return _res;
}

static PyObject *glut_glutDeviceGet(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    int _rv;
    long type;
    if (!PyArg_ParseTuple(_args, "l",
			  &type))
	return NULL;
    _rv = glutDeviceGet(type);
    _res = Py_BuildValue("i",
			 _rv);
    return _res;
}

static PyObject *glut_glutExtensionSupported(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    int _rv;
    char *name;
    if (!PyArg_ParseTuple(_args, "s",
			  &name))
	return NULL;
    _rv = glutExtensionSupported(name);
    _res = Py_BuildValue("i",
			 _rv);
    return _res;
}

#if (GLUT_API_VERSION >= 3)
static PyObject *glut_glutGetModifiers(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    int _rv;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    _rv = glutGetModifiers();
    _res = Py_BuildValue("i",
			 _rv);
    return _res;
}

static PyObject *glut_glutLayerGet(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    int _rv;
    long type;
    if (!PyArg_ParseTuple(_args, "l",
			  &type))
	return NULL;
    _rv = glutLayerGet(type);
    _res = Py_BuildValue("i",
			 _rv);
    return _res;
}

#endif

static PyObject *glut_glutBitmapCharacter(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    char *font_name;
    void *font;
    int character;
    if (!PyArg_ParseTuple(_args, "si",
			  &font_name,
			  &character))
	return NULL;
    font = map_lookup(glutFonts, font_name);
    glutBitmapCharacter(font,
			character);
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutBitmapWidth(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    int _rv;
    char *font_name;
    void *font;
    int character;
    if (!PyArg_ParseTuple(_args, "si",
			  &font_name,
			  &character))
	return NULL;
    font = map_lookup(glutFonts, font_name);
    _rv = glutBitmapWidth(font,
			  character);
    _res = Py_BuildValue("i",
			 _rv);
    return _res;
}

static PyObject *glut_glutStrokeCharacter(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    char *font_name;
    void *font;
    int character;
    if (!PyArg_ParseTuple(_args, "si",
			  &font_name,
			  &character))
	return NULL;
    font = map_lookup(glutFonts, font_name);
    glutStrokeCharacter(font,
			character);
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutStrokeWidth(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    int _rv;
    char *font_name;
    void *font;
    int character;
    if (!PyArg_ParseTuple(_args, "si",
			  &font_name,
			  &character))
	return NULL;
    font = map_lookup(glutFonts, font_name);
    _rv = glutStrokeWidth(font,
			  character);
    _res = Py_BuildValue("i",
			 _rv);
    return _res;
}

static PyObject *glut_glutWireSphere(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    double radius;
    int slices;
    int stacks;
    if (!PyArg_ParseTuple(_args, "dii",
			  &radius,
			  &slices,
			  &stacks))
	return NULL;
    glutWireSphere(radius,
		   slices,
		   stacks);
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutSolidSphere(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    double radius;
    int slices;
    int stacks;
    if (!PyArg_ParseTuple(_args, "dii",
			  &radius,
			  &slices,
			  &stacks))
	return NULL;
    glutSolidSphere(radius,
		    slices,
		    stacks);
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutWireCone(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    double base;
    double height;
    int slices;
    int stacks;
    if (!PyArg_ParseTuple(_args, "ddii",
			  &base,
			  &height,
			  &slices,
			  &stacks))
	return NULL;
    glutWireCone(base,
		 height,
		 slices,
		 stacks);
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutSolidCone(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    double base;
    double height;
    int slices;
    int stacks;
    if (!PyArg_ParseTuple(_args, "ddii",
			  &base,
			  &height,
			  &slices,
			  &stacks))
	return NULL;
    glutSolidCone(base,
		  height,
		  slices,
		  stacks);
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutWireCube(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    double size;
    if (!PyArg_ParseTuple(_args, "d",
			  &size))
	return NULL;
    glutWireCube(size);
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutSolidCube(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    double size;
    if (!PyArg_ParseTuple(_args, "d",
			  &size))
	return NULL;
    glutSolidCube(size);
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

/* These functions are just like the solidCube and wireCube, but for
   boxes with any dimension */

/* this is taken from glut_shapes.c */
static void drawBox(GLfloat x0, GLfloat x1, GLfloat y0, GLfloat y1,
		    GLfloat z0, GLfloat z1, GLenum type)
{
    static GLfloat n[6][3] =
    {
	{-1.0, 0.0, 0.0},
	{0.0, 1.0, 0.0},
	{1.0, 0.0, 0.0},
	{0.0, -1.0, 0.0},
	{0.0, 0.0, 1.0},
	{0.0, 0.0, -1.0}
    };
    static GLint faces[6][4] =
    {
	{0, 1, 2, 3},
	{3, 2, 6, 7},
	{7, 6, 5, 4},
	{4, 5, 1, 0},
	{5, 6, 2, 1},
	{7, 4, 0, 3}
    };
    GLfloat v[8][3], tmp;
    GLint i;

    if (x0 > x1) {
	tmp = x0;
	x0 = x1;
	x1 = tmp;
    }
    if (y0 > y1) {
	tmp = y0;
	y0 = y1;
	y1 = tmp;
    }
    if (z0 > z1) {
	tmp = z0;
	z0 = z1;
	z1 = tmp;
    }
    v[0][0] = v[1][0] = v[2][0] = v[3][0] = x0;
    v[4][0] = v[5][0] = v[6][0] = v[7][0] = x1;
    v[0][1] = v[1][1] = v[4][1] = v[5][1] = y0;
    v[2][1] = v[3][1] = v[6][1] = v[7][1] = y1;
    v[0][2] = v[3][2] = v[4][2] = v[7][2] = z0;
    v[1][2] = v[2][2] = v[5][2] = v[6][2] = z1;

    for (i = 0; i < 6; i++) {
	glBegin(type);
	glNormal3fv(&n[i][0]);
	glVertex3fv(&v[faces[i][0]][0]);
	glVertex3fv(&v[faces[i][1]][0]);
	glVertex3fv(&v[faces[i][2]][0]);
	glVertex3fv(&v[faces[i][3]][0]);
	glEnd();
    }
}

static void glutWireBox(GLdouble x0, GLdouble x1, GLdouble y0, GLdouble y1,
			GLdouble z0, GLdouble z1)
{
    drawBox(x0, x1, y0, y1, z0, z1, GL_LINE_LOOP);
}

static void glutSolidBox(GLdouble x0, GLdouble y0, GLdouble z0, GLdouble x1,
			 GLdouble y1, GLdouble z1)
{
    drawBox(x0, x1, y0, y1, z0, z1, GL_QUADS);
}
static PyObject *glut_glutWireBox(PyObject * self, PyObject * args)
{
    PyObject *res = NULL;
    double x1, y1, z1, x2, y2, z2;
    if (!PyArg_ParseTuple(args, "dddddd", &x1, &y1, &z1, &x2, &y2, &z2))
	return NULL;
    glutWireBox(x1, x2, y1, y2, z1, z2);
    Py_INCREF(Py_None);
    res = Py_None;
    return res;
}

static PyObject *glut_glutSolidBox(PyObject * self, PyObject * args)
{
    PyObject *res = NULL;
    double x1, y1, z1, x2, y2, z2;
    if (!PyArg_ParseTuple(args, "dddddd", &x1, &y1, &z1, &x2, &y2, &z2))
	return NULL;
    glutSolidBox(x1, x2, y1, y2, z1, z2);
    Py_INCREF(Py_None);
    res = Py_None;
    return res;
}

static PyObject *glut_glutWireTorus(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    double innerRadius;
    double outerRadius;
    int sides;
    int rings;
    if (!PyArg_ParseTuple(_args, "ddii",
			  &innerRadius,
			  &outerRadius,
			  &sides,
			  &rings))
	return NULL;
    glutWireTorus(innerRadius,
		  outerRadius,
		  sides,
		  rings);
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutSolidTorus(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    double innerRadius;
    double outerRadius;
    int sides;
    int rings;
    if (!PyArg_ParseTuple(_args, "ddii",
			  &innerRadius,
			  &outerRadius,
			  &sides,
			  &rings))
	return NULL;
    glutSolidTorus(innerRadius,
		   outerRadius,
		   sides,
		   rings);
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutWireDodecahedron(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutWireDodecahedron();
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutSolidDodecahedron(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutSolidDodecahedron();
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutWireTeapot(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    double size;
    if (!PyArg_ParseTuple(_args, "d",
			  &size))
	return NULL;
    glutWireTeapot(size);
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutSolidTeapot(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    double size;
    if (!PyArg_ParseTuple(_args, "d",
			  &size))
	return NULL;
    glutSolidTeapot(size);
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutWireOctahedron(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutWireOctahedron();
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutSolidOctahedron(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutSolidOctahedron();
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutWireTetrahedron(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutWireTetrahedron();
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutSolidTetrahedron(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutSolidTetrahedron();
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutWireIcosahedron(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutWireIcosahedron();
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

static PyObject *glut_glutSolidIcosahedron(PyObject * _self, PyObject * _args)
{
    PyObject *_res = NULL;
    if (!PyArg_ParseTuple(_args, ""))
	return NULL;
    glutSolidIcosahedron();
    Py_INCREF(Py_None);
    _res = Py_None;
    return _res;
}

#ifdef SGISTEREO

#include <X11/Xatom.h>
#include <X11/extensions/SGIStereo.h>

static PyObject *stereo_Error;

/* Standard screen dimensions */
#define XMAXSCREEN      1280
#define YMAXSCREEN      1024

#define YSTEREO         491	/* Subfield height in pixels */
#define YOFFSET_LEFT    532	/* YSTEREO + YBLANK */

static char stereo_start__doc__[] = "";
static PyObject *
 stereo_start(PyObject * self, PyObject * args)
{
    int event, error;

    if (!PyArg_ParseTuple(args, ""))
	return NULL;
    if (!XSGIStereoQueryExtension(__glutDisplay, &event, &error)) {
	fprintf(stderr, "Stereo not supported on this display!\n");
	exit(0);
    }
    if (XSGIQueryStereoMode(__glutDisplay, __glutCurrentWindow->win) < 0) {
	fprintf(stderr, "Stereo not supported on this window!\n");
	exit(0);
    }
    if (system("/usr/gfx/setmon -n STR_BOT") != 0) {
	fprintf(stderr, "setmon attempt failed!\n");
	system("/usr/gfx/setmon -n 72hz");
	exit(0);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static char stereo_stop__doc__[] = "";

static PyObject *
 stereo_stop(PyObject * self, PyObject * args)
{

    if (!PyArg_ParseTuple(args, ""))
	return NULL;
    system("/usr/gfx/setmon -n 72hz");

    Py_INCREF(Py_None);
    return Py_None;
}

static char stereo_left__doc__[] = "";
static PyObject *
 stereo_left(PyObject * self, PyObject * args)
{
    if (!PyArg_ParseTuple(args, ""))
	return NULL;
    XSGISetStereoBuffer(__glutDisplay, __glutCurrentWindow->win, STEREO_BUFFER_LEFT);
    glViewport(0, 0, XMAXSCREEN, YSTEREO);

    XSync(__glutDisplay, False);

    Py_INCREF(Py_None);
    return Py_None;
}

static char stereo_right__doc__[] = "";
static PyObject *
 stereo_right(PyObject * self, PyObject * args)
{
    if (!PyArg_ParseTuple(args, ""))
	return NULL;

    XSGISetStereoBuffer(__glutDisplay, __glutCurrentWindow->win, STEREO_BUFFER_RIGHT);
    glViewport(0, 0, XMAXSCREEN, YSTEREO);

    XSync(__glutDisplay, False);

    Py_INCREF(Py_None);
    return Py_None;
}
#endif

static PyMethodDef glut_methods[] =
{
    {"glutSetCreateMenuCallback", (PyCFunction) glut_SetglutCreateMenuCallback, 1},
    {"glutSetDisplayFuncCallback", (PyCFunction) glut_SetglutDisplayFuncCallback, 1},
    {"glutSetReshapeFuncCallback", (PyCFunction) glut_SetglutReshapeFuncCallback, 1},
    {"glutSetKeyboardFuncCallback", (PyCFunction) glut_SetglutKeyboardFuncCallback, 1},
#if (GLUT_API_VERSION >= 4)
    {"glutSetKeyboardUpFuncCallback", (PyCFunction) glut_SetglutKeyboardUpFuncCallback, 1},
#endif
    {"glutSetMouseFuncCallback", (PyCFunction) glut_SetglutMouseFuncCallback, 1},
    {"glutSetMotionFuncCallback", (PyCFunction) glut_SetglutMotionFuncCallback, 1},
    {"glutSetPassiveMotionFuncCallback", (PyCFunction) glut_SetglutPassiveMotionFuncCallback, 1},
    {"glutSetEntryFuncCallback", (PyCFunction) glut_SetglutEntryFuncCallback, 1},
    {"glutSetVisibilityFuncCallback", (PyCFunction) glut_SetglutVisibilityFuncCallback, 1},
    {"glutSetIdleFuncCallback", (PyCFunction) glut_SetglutIdleFuncCallback, 1},
    {"glutSetTimerFuncCallback", (PyCFunction) glut_SetglutTimerFuncCallback, 1},
    {"glutSetMenuStateFuncCallback", (PyCFunction) glut_SetglutMenuStateFuncCallback, 1},
    {"glutSetSpecialFuncCallback", (PyCFunction) glut_SetglutSpecialFuncCallback, 1},
#if (GLUT_API_VERSION >= 4)
    {"glutSetSpecialUpFuncCallback", (PyCFunction) glut_SetglutSpecialUpFuncCallback, 1},
#endif
    {"glutSetSpaceballMotionFuncCallback", (PyCFunction) glut_SetglutSpaceballMotionFuncCallback, 1},
    {"glutSetSpaceballRotateFuncCallback", (PyCFunction) glut_SetglutSpaceballRotateFuncCallback, 1},
    {"glutSetSpaceballButtonFuncCallback", (PyCFunction) glut_SetglutSpaceballButtonFuncCallback, 1},
    {"glutSetButtonBoxFuncCallback", (PyCFunction) glut_SetglutButtonBoxFuncCallback, 1},
    {"glutSetDialsFuncCallback", (PyCFunction) glut_SetglutDialsFuncCallback, 1},
    {"glutSetTabletMotionFuncCallback", (PyCFunction) glut_SetglutTabletMotionFuncCallback, 1},
    {"glutSetTabletButtonFuncCallback", (PyCFunction) glut_SetglutTabletButtonFuncCallback, 1},
    {"glutSetMenuStatusFuncCallback", (PyCFunction) glut_SetglutMenuStatusFuncCallback, 1},
    {"glutSetOverlayDisplayFuncCallback", (PyCFunction) glut_SetglutOverlayDisplayFuncCallback, 1},
    {"glutInit", (PyCFunction) glut_glutInit, 1,
     "(char* argv) -> (int argcp)"},
    {"glutInitDisplayMode", (PyCFunction) glut_glutInitDisplayMode, 1,
     "(int mode) -> None"},
    {"glutInitWindowPosition", (PyCFunction) glut_glutInitWindowPosition, 1,
     "(int x, int y) -> None"},
    {"glutInitWindowSize", (PyCFunction) glut_glutInitWindowSize, 1,
     "(int width, int height) -> None"},
    {"glutMainLoop", (PyCFunction) glut_glutMainLoop, 1,
     "() -> None"},
    {"glutCreateWindow", (PyCFunction) glut_glutCreateWindow, 1,
     "(char* title) -> (int _rv)"},
    {"glutCreateSubWindow", (PyCFunction) glut_glutCreateSubWindow, 1,
     "(int win, int x, int y, int width, int height) -> (int _rv)"},
    {"glutDestroyWindow", (PyCFunction) glut_glutDestroyWindow, 1,
     "(int win) -> None"},
    {"glutPostRedisplay", (PyCFunction) glut_glutPostRedisplay, 1,
     "() -> None"},
    {"glutSwapBuffers", (PyCFunction) glut_glutSwapBuffers, 1,
     "() -> None"},
    {"glutGetWindow", (PyCFunction) glut_glutGetWindow, 1,
     "() -> (int _rv)"},
    {"glutSetWindow", (PyCFunction) glut_glutSetWindow, 1,
     "(int win) -> None"},
    {"glutSetWindowTitle", (PyCFunction) glut_glutSetWindowTitle, 1,
     "(char* title) -> None"},
    {"glutSetIconTitle", (PyCFunction) glut_glutSetIconTitle, 1,
     "(char* title) -> None"},
    {"glutPositionWindow", (PyCFunction) glut_glutPositionWindow, 1,
     "(int x, int y) -> None"},
    {"glutReshapeWindow", (PyCFunction) glut_glutReshapeWindow, 1,
     "(int width, int height) -> None"},
    {"glutPopWindow", (PyCFunction) glut_glutPopWindow, 1,
     "() -> None"},
    {"glutPushWindow", (PyCFunction) glut_glutPushWindow, 1,
     "() -> None"},
    {"glutIconifyWindow", (PyCFunction) glut_glutIconifyWindow, 1,
     "() -> None"},
    {"glutShowWindow", (PyCFunction) glut_glutShowWindow, 1,
     "() -> None"},
    {"glutHideWindow", (PyCFunction) glut_glutHideWindow, 1,
     "() -> None"},
    {"glutFullScreen", (PyCFunction) glut_glutFullScreen, 1,
     "() -> None"},
#if (GLUT_API_VERSION >= 4)
    {"glutGameModeString", (PyCFunction) glut_glutGameModeString, 1,
     "(char* s) -> None"},
    {"glutEnterGameMode", (PyCFunction) glut_glutEnterGameMode, 1,
     "() -> None"},
    {"glutLeaveGameMode", (PyCFunction) glut_glutLeaveGameMode, 1,
     "() -> None"},
    {"glutGameModeGet", (PyCFunction) glut_glutGameModeGet, 1,
     "(long type) -> (int _rv)"},
#endif
    {"glutSetCursor", (PyCFunction) glut_glutSetCursor, 1,
     "(int cursor) -> None"},
#if (GLUT_API_VERSION >= 4)
    {"glutIgnoreKeyRepeat", (PyCFunction) glut_glutIgnoreKeyRepeat, 1,
     "(int value) -> None"},
#endif
    {"glutEstablishOverlay", (PyCFunction) glut_glutEstablishOverlay, 1,
     "() -> None"},
    {"glutRemoveOverlay", (PyCFunction) glut_glutRemoveOverlay, 1,
     "() -> None"},
    {"glutUseLayer", (PyCFunction) glut_glutUseLayer, 1,
     "(long layer) -> None"},
    {"glutPostOverlayRedisplay", (PyCFunction) glut_glutPostOverlayRedisplay, 1,
     "() -> None"},
    {"glutShowOverlay", (PyCFunction) glut_glutShowOverlay, 1,
     "() -> None"},
    {"glutHideOverlay", (PyCFunction) glut_glutHideOverlay, 1,
     "() -> None"},
    {"glutCreateMenu", (PyCFunction) glut_glutCreateMenu, 1,
     "(staticGlobal glutCreateMenuCallback) -> (int _rv)"},
    {"glutDestroyMenu", (PyCFunction) glut_glutDestroyMenu, 1,
     "(int menu) -> None"},
    {"glutGetMenu", (PyCFunction) glut_glutGetMenu, 1,
     "() -> (int _rv)"},
    {"glutSetMenu", (PyCFunction) glut_glutSetMenu, 1,
     "(int menu) -> None"},
    {"glutAddMenuEntry", (PyCFunction) glut_glutAddMenuEntry, 1,
     "(char* label, int value) -> None"},
    {"glutAddSubMenu", (PyCFunction) glut_glutAddSubMenu, 1,
     "(char* label, int submenu) -> None"},
    {"glutChangeToMenuEntry", (PyCFunction) glut_glutChangeToMenuEntry, 1,
     "(int item, char* label, int value) -> None"},
    {"glutChangeToSubMenu", (PyCFunction) glut_glutChangeToSubMenu, 1,
     "(int item, char* label, int submenu) -> None"},
    {"glutRemoveMenuItem", (PyCFunction) glut_glutRemoveMenuItem, 1,
     "(int item) -> None"},
    {"glutAttachMenu", (PyCFunction) glut_glutAttachMenu, 1,
     "(int button) -> None"},
    {"glutDetachMenu", (PyCFunction) glut_glutDetachMenu, 1,
     "(int button) -> None"},
    {"glutDisplayFunc", (PyCFunction) glut_glutDisplayFunc, 1,
     "(staticGlobal glutDisplayFuncCallback) -> None"},
    {"glutReshapeFunc", (PyCFunction) glut_glutReshapeFunc, 1,
     "(staticGlobal glutReshapeFuncCallback) -> None"},
    {"glutKeyboardFunc", (PyCFunction) glut_glutKeyboardFunc, 1,
     "(staticGlobal glutKeyboardFuncCallback) -> None"},
#if (GLUT_API_VERSION >= 4)
    {"glutKeyboardUpFunc", (PyCFunction) glut_glutKeyboardUpFunc, 1,
     "(staticGlobal glutKeyboardUpFuncCallback) -> None"},
#endif
    {"glutMouseFunc", (PyCFunction) glut_glutMouseFunc, 1,
     "(staticGlobal glutMouseFuncCallback) -> None"},
    {"glutMotionFunc", (PyCFunction) glut_glutMotionFunc, 1,
     "(staticGlobal glutMotionFuncCallback) -> None"},
    {"glutPassiveMotionFunc", (PyCFunction) glut_glutPassiveMotionFunc, 1,
     "(staticGlobal glutPassiveMotionFuncCallback) -> None"},
    {"glutEntryFunc", (PyCFunction) glut_glutEntryFunc, 1,
     "(staticGlobal glutEntryFuncCallback) -> None"},
    {"glutVisibilityFunc", (PyCFunction) glut_glutVisibilityFunc, 1,
     "(staticGlobal glutVisibilityFuncCallback) -> None"},
    {"glutIdleFunc", (PyCFunction) glut_glutIdleFunc, 1,
     "(staticGlobal glutIdleFuncCallback) -> None"},
    {"glutTimerFunc", (PyCFunction) glut_glutTimerFunc, 1,
  "(int millis, staticGlobal glutTimerFuncCallback, int value) -> None"},
    {"glutMenuStateFunc", (PyCFunction) glut_glutMenuStateFunc, 1,
     "(staticGlobal glutMenuStateFuncCallback) -> None"},
    {"glutSpecialFunc", (PyCFunction) glut_glutSpecialFunc, 1,
     "(staticGlobal glutSpecialFuncCallback) -> None"},
#if (GLUT_API_VERSION >= 4)
    {"glutSpecialUpFunc", (PyCFunction) glut_glutSpecialUpFunc, 1,
     "(staticGlobal glutSpecialFuncUpCallback) -> None"},
#endif
{"glutSpaceballMotionFunc", (PyCFunction) glut_glutSpaceballMotionFunc, 1,
 "(staticGlobal glutSpaceballMotionFuncCallback) -> None"},
{"glutSpaceballRotateFunc", (PyCFunction) glut_glutSpaceballRotateFunc, 1,
 "(staticGlobal glutSpaceballRotateFuncCallback) -> None"},
{"glutSpaceballButtonFunc", (PyCFunction) glut_glutSpaceballButtonFunc, 1,
 "(staticGlobal glutSpaceballButtonFuncCallback) -> None"},
    {"glutButtonBoxFunc", (PyCFunction) glut_glutButtonBoxFunc, 1,
     "(staticGlobal glutButtonBoxFuncCallback) -> None"},
    {"glutDialsFunc", (PyCFunction) glut_glutDialsFunc, 1,
     "(staticGlobal glutDialsFuncCallback) -> None"},
    {"glutTabletMotionFunc", (PyCFunction) glut_glutTabletMotionFunc, 1,
     "(staticGlobal glutTabletMotionFuncCallback) -> None"},
    {"glutTabletButtonFunc", (PyCFunction) glut_glutTabletButtonFunc, 1,
     "(staticGlobal glutTabletButtonFuncCallback) -> None"},
    {"glutMenuStatusFunc", (PyCFunction) glut_glutMenuStatusFunc, 1,
     "(staticGlobal glutMenuStatusFuncCallback) -> None"},
 {"glutOverlayDisplayFunc", (PyCFunction) glut_glutOverlayDisplayFunc, 1,
  "(staticGlobal glutOverlayDisplayFuncCallback) -> None"},
    {"glutSetColor", (PyCFunction) glut_glutSetColor, 1,
     "(int x, float red, float green, float blue) -> None"},
    {"glutGetColor", (PyCFunction) glut_glutGetColor, 1,
     "(int ndx, int component) -> (float _rv)"},
    {"glutCopyColormap", (PyCFunction) glut_glutCopyColormap, 1,
     "(int win) -> None"},
    {"glutGet", (PyCFunction) glut_glutGet, 1,
     "(long type) -> (int _rv)"},
    {"glutDeviceGet", (PyCFunction) glut_glutDeviceGet, 1,
     "(long type) -> (int _rv)"},
 {"glutExtensionSupported", (PyCFunction) glut_glutExtensionSupported, 1,
  "(char* name) -> (int _rv)"},
#if (GLUT_API_VERSION >= 3)
    {"glutGetModifiers", (PyCFunction) glut_glutGetModifiers, 1,
     "() -> (int _rv)"},
    {"glutLayerGet", (PyCFunction) glut_glutLayerGet, 1,
     "(long type) -> (int _rv)"},
#endif
    {"glutBitmapCharacter", (PyCFunction) glut_glutBitmapCharacter, 1,
     "(stringMap font, int character) -> None"},
    {"glutBitmapWidth", (PyCFunction) glut_glutBitmapWidth, 1,
     "(stringMap font, int character) -> (int _rv)"},
    {"glutStrokeCharacter", (PyCFunction) glut_glutStrokeCharacter, 1,
     "(stringMap font, int character) -> None"},
    {"glutStrokeWidth", (PyCFunction) glut_glutStrokeWidth, 1,
     "(stringMap font, int character) -> (int _rv)"},
    {"glutWireSphere", (PyCFunction) glut_glutWireSphere, 1,
     "(double radius, int slices, int stacks) -> None"},
    {"glutSolidSphere", (PyCFunction) glut_glutSolidSphere, 1,
     "(double radius, int slices, int stacks) -> None"},
    {"glutWireCone", (PyCFunction) glut_glutWireCone, 1,
     "(double base, double height, int slices, int stacks) -> None"},
    {"glutSolidCone", (PyCFunction) glut_glutSolidCone, 1,
     "(double base, double height, int slices, int stacks) -> None"},
    {"glutWireCube", (PyCFunction) glut_glutWireCube, 1,
     "(double size) -> None"},
    {"glutSolidCube", (PyCFunction) glut_glutSolidCube, 1,
     "(double size) -> None"},
    {"glutWireBox", glut_glutWireBox, 1,
     "(double x1, double y1, double z1, double x2, double y2, double z2) -> None"},
    {"glutSolidBox", glut_glutSolidBox, 1,
     "(double x1, double y1, double z1, double x2, double y2, double z2) -> None"},
    {"glutWireTorus", (PyCFunction) glut_glutWireTorus, 1,
"(double innerRadius, double outerRadius, int sides, int rings) -> None"},
    {"glutSolidTorus", (PyCFunction) glut_glutSolidTorus, 1,
"(double innerRadius, double outerRadius, int sides, int rings) -> None"},
    {"glutWireDodecahedron", (PyCFunction) glut_glutWireDodecahedron, 1,
     "() -> None"},
    {"glutSolidDodecahedron", (PyCFunction) glut_glutSolidDodecahedron, 1,
     "() -> None"},
    {"glutWireTeapot", (PyCFunction) glut_glutWireTeapot, 1,
     "(double size) -> None"},
    {"glutSolidTeapot", (PyCFunction) glut_glutSolidTeapot, 1,
     "(double size) -> None"},
    {"glutWireOctahedron", (PyCFunction) glut_glutWireOctahedron, 1,
     "() -> None"},
    {"glutSolidOctahedron", (PyCFunction) glut_glutSolidOctahedron, 1,
     "() -> None"},
    {"glutWireTetrahedron", (PyCFunction) glut_glutWireTetrahedron, 1,
     "() -> None"},
    {"glutSolidTetrahedron", (PyCFunction) glut_glutSolidTetrahedron, 1,
     "() -> None"},
    {"glutWireIcosahedron", (PyCFunction) glut_glutWireIcosahedron, 1,
     "() -> None"},
    {"glutSolidIcosahedron", (PyCFunction) glut_glutSolidIcosahedron, 1,
     "() -> None"},
#ifdef SGISTEREO
    {"glutstart_stereo", stereo_start, 1, stereo_start__doc__},
    {"glutstop_stereo", stereo_stop, 1, stereo_stop__doc__},
    {"glutleft", stereo_left, 1, stereo_left__doc__},
    {"glutright", stereo_right, 1, stereo_right__doc__},
#endif
    {NULL, NULL, 0}
};

DL_EXPORT(void) init_glut(void);
DL_EXPORT(void) init_glut(void)
{
    PyObject *m;
    PyObject *d;

    m = Py_InitModule("_glut", glut_methods);
    d = PyModule_GetDict(m);
    glut_Error = PyString_FromString("_glut.Error");
    if (glut_Error == NULL ||
	PyDict_SetItemString(d, "Error", glut_Error) != 0)
	Py_FatalError("can't initialize glut.Error");
}
#else

/* NO GLUT available */

static PyMethodDef glut_methods[] =
{
    {NULL, NULL, 0}
};

DL_EXPORT(void)
init_glut(void)
{
    PyObject *m;
    PyObject *d;

    m = Py_InitModule("_glut", glut_methods);
    d = PyModule_GetDict(m);
}

#endif

#ifdef WIN32
/* for distutils compatibility on WIN32 */
DL_EXPORT(void) init_glutmodule(void) {init_glut();}
#endif
