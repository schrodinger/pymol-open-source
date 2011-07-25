#ifdef _PYMOL_WIN32
#include"os_predef.h"
#endif
/**
 * $Id$
 *
 * OpenGL Module for Python
 * 
 * Version: 0.9.PyMOL (modified for usage with PyMOL)
 *
 * Date:    1998/01/15
 *
 * Authors: Tom Schwaller     <tom.schwaller@linux-magazin.de>
 *          Jim Hugunin       <hugunin@python.org>
 *          David Ascher      <da@skivs.ski.org>
 *          Michel Sanner     <sanner@scripps.edu>
 * 
 * Notes:
 *
 *   - the file glconst.py is also needed, as it defines all kinds of
 *     useful constants
 *     
 *   - If you do not have Numeric python installed, undefine NUMERIC
 *   - We have included the glut shapes for convenience. If you do not 
 *     want to use them undefine the GLUT constant, otherwise add 
 *     glut_shapes.c in the setup file like:
 *     opengl openglmodule.c glut_shapes.c -lGL -lGLU
 *
 *   - thanks to Mike Hartshorn <mjh@yorvic.york.ac.uk> for quite 
 *     some bug reports.
 *
 *   - Some new stuff: Hope to have more from users of this module
 *
 *     1996/08/30: py_gl_TranslateScene, py_gl_RotateScene
 *     1996/09/01: py_gl_DistFromLine (from Mike)
 * 
 *     1997/03/02: py_gl_SaveTiff from Lothar Birk <lb@ocean.fb12.TU-Berlin.DE>
 *                 use #define LIBTIFF 1 or -DLIBTIFF on the comand line
 *                 and -ltiff as axtra library
 *
 *     1997/05/12: py_gl_SavePPM fixed by Tom Schwaller 
 *     1997/05/12: py_gl_SaveEPS added by Tom Schwaller 
 *                 code from Togl1.3, was contributed by
 *                 Miguel A. De Riera Pasenau <miguel@DALILA.UPC.ES>
 *
 *     1998/01/16: moved all non-OpenGL spec stuff into openglutil.c
 *                 [DA & MS]
 *
 * 
***/

#ifdef WIN32
#define MS_WIN32
#endif

#include "Python.h"
#ifdef MS_WIN32
#include <windows.h>
#endif

#ifndef _PYMOL_OSX
#include <GL/gl.h>
/* is there a way to check and see if we have glu? */
#else
#include <gl.h>
#endif

#include <math.h>

#ifdef NUMERIC
static int typecode2gltype[] =
    { -1, GL_UNSIGNED_BYTE, GL_BYTE, GL_SHORT, GL_INT, -1, GL_FLOAT, -1,
    -1, -1, -1, -1
};
#endif

static PyObject *py_gl_Error;

static PyObject *ErrorReturn(char *message)
{
    PyErr_SetString(py_gl_Error, message);
    return NULL;
}

#define ASSERT(E,M) if(!(E)) return ErrorReturn(M)

#define TRY(E) if(! (E)) return NULL

#ifdef NUMERIC

#if   defined(HAVE_ARRAYOBJECT_H)
#include "arrayobject.h"
#elif defined(HAVE_EXTENSIONS_ARRAYOBJECT_H)
#include "Extensions/arrayobject.h"
#elif defined(HAVE_NUMERIC_ARRAYOBJECT_H)
#include "numpy/oldnumeric.h"
#elif defined(HAVE_NUMERICAL_ARRAYOBJECT_H)
#include "numerical/arrayobject.h"
#else
/*#error "Don't know where to find file 'arrayobject.h'" */
#include "arrayobject.h"
#endif

#define PyArray_ClearMemory(op, pitems) Py_DECREF((op))

static int PyArray_AsUByteArray(PyObject ** op, GLubyte ** pitems, int *pn)
{
    PyArrayObject *mp;
    mp =
	(PyArrayObject *) PyArray_ContiguousFromObject(*op, PyArray_UBYTE,
						       0, 1);
    if (!mp)
	return 0;
    *pitems = (GLubyte *) (mp->data);
    *pn = mp->dimensions[0];
    *op = (PyObject *) mp;
    return 1;
}

static int PyArray_AsDoubleArray(PyObject ** op, GLdouble ** pitems, int *pn)
{
    PyArrayObject *mp;
    mp =
	(PyArrayObject *) PyArray_ContiguousFromObject(*op,
						       PyArray_DOUBLE, 0, 1);
    if (!mp)
	return 0;
    *pitems = (GLdouble *) (mp->data);
    *pn = mp->dimensions[0];
    *op = (PyObject *) mp;
    return 1;
}

static int PyArray_AsFloatArray(PyObject ** op, GLfloat ** pitems, int *pn)
{
    PyArrayObject *mp;
    mp =
	(PyArrayObject *) PyArray_ContiguousFromObject(*op, PyArray_FLOAT,
						       0, 1);
    if (!mp)
	return 0;
    *pitems = (GLfloat *) (mp->data);
    *pn = mp->dimensions[0];
    *op = (PyObject *) mp;
    return 1;
}

static int PyArray_AsIntArray(PyObject ** op, GLint ** pitems, int *pn)
{
    PyArrayObject *mp;
    mp = (PyArrayObject *) PyArray_ContiguousFromObject(*op, PyArray_INT, 0, 1);
    if (!mp)
	return 0;
    *pitems = (GLint *) (mp->data);
    *pn = mp->dimensions[0];
    *op = (PyObject *) mp;
    return 1;
}

static int PyArray_AsShortArray(PyObject ** op, GLshort ** pitems, int *pn)
{
    PyArrayObject *mp;
    mp =
	(PyArrayObject *) PyArray_ContiguousFromObject(*op, PyArray_SHORT,
						       0, 1);
    if (!mp)
	return 0;
    *pitems = (GLshort *) (mp->data);
    *pn = mp->dimensions[0];
    *op = (PyObject *) mp;
    return 1;
}

#else

#include "abstract.h"

#define PyArray_ClearMemory(op, pitems) PyMem_DEL(pitems) /* is pitems PyMem or PyObject? */

static int PyArray_AsDoubleArray(PyObject ** op, GLdouble ** pitems, int *pn)
{
    GLdouble *items;
    PyObject *item;
    int n, i;
    if (!PySequence_Check(*op))
	return 0;
    n = PySequence_Length(*op);
    items = PyMem_NEW(GLdouble, n);
    if (items == NULL) {
	PyErr_NoMemory();
	return 0;
    }
    for (i = 0; i < n; i++) {
	if ((item = PySequence_GetItem(*op, i))) {
	    items[i] = PyFloat_AsDouble(item);
	    Py_DECREF(item);
	}
	if (PyErr_Occurred())
	    return 0;
    }
    *pitems = items;
    *pn = n;
    return 1;
}

static int PyArray_AsFloatArray(PyObject ** op, GLfloat ** pitems, int *pn)
{
    GLfloat *items;
    PyObject *item;
    int n, i;
    if (!PySequence_Check(*op))
	return 0;
    n = PySequence_Length(*op);
    items = PyMem_NEW(GLfloat, n);
    if (items == NULL) {
	PyErr_NoMemory();
	return 0;
    }
    for (i = 0; i < n; i++) {
	if ((item = PySequence_GetItem(*op, i))) {
	    items[i] = (float)PyFloat_AsDouble(item);
	    Py_DECREF(item);
	}
	if (PyErr_Occurred())
	    return 0;
    }
    *pitems = items;
    *pn = n;
    return 1;
}

static int PyArray_AsIntArray(PyObject ** op, GLint ** pitems, int *pn)
{
    GLint *items;
    PyObject *item;
    int n, i;
    if (!PySequence_Check(*op))
	return 0;
    n = PySequence_Length(*op);
    items = PyMem_NEW(GLint, n);
    if (items == NULL) {
	PyErr_NoMemory();
	return 0;
    }
    for (i = 0; i < n; i++) {
	if ((item = PySequence_GetItem(*op, i))) {
	    items[i] = PyInt_AsLong(item);
	    Py_DECREF(item);
	}
	if (PyErr_Occurred())
	    return 0;
    }
    *pitems = items;
    *pn = n;
    return 1;
}

static int PyArray_AsShortArray(PyObject ** op, GLshort ** pitems, int *pn)
{
    GLshort *items;
    PyObject *item;
    int n, i;
    if (!PySequence_Check(*op))
	return 0;
    n = PySequence_Length(*op);
    items = PyMem_NEW(GLshort, n);
    if (items == NULL) {
	PyErr_NoMemory();
	return 0;
    }
    for (i = 0; i < n; i++) {
	if ((item = PySequence_GetItem(*op, i))) {
	    items[i] = (GLshort)PyInt_AsLong(item);
	    Py_DECREF(item);
	}
	if (PyErr_Occurred())
	    return 0;
    }
    *pitems = items;
    *pn = n;
    return 1;

}

#endif				/* Not NUMERIC */

/* 
   ########################################################################
 */

static PyObject *py_gl_Accum(PyObject * self, PyObject * args)
{
    GLenum arg1;
    GLfloat arg2;
    TRY(PyArg_ParseTuple(args, "if", &arg1, &arg2));
    glAccum(arg1, arg2);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_AlphaFunc(PyObject * self, PyObject * args)
{
    GLenum arg1;
    GLclampf arg2;
    TRY(PyArg_ParseTuple(args, "if", &arg1, &arg2));
    glAlphaFunc(arg1, arg2);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Begin(PyObject * self, PyObject * args)
{
    GLenum arg1;
    TRY(PyArg_ParseTuple(args, "i", &arg1));
    glBegin(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_BlendFunc(PyObject * self, PyObject * args)
{
    GLenum arg1, arg2;
    TRY(PyArg_ParseTuple(args, "ii", &arg1, &arg2));
    glBlendFunc(arg1, arg2);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_CallList(PyObject * self, PyObject * args)
{
    GLuint arg1;
    TRY(PyArg_ParseTuple(args, "i", &arg1));
    glCallList(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Clear(PyObject * self, PyObject * args)
{
    GLbitfield arg1;
    TRY(PyArg_ParseTuple(args, "i", &arg1));
    glClear(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_ClearAccum(PyObject * self, PyObject * args)
{
    GLfloat arg1, arg2, arg3, arg4;
    TRY(PyArg_ParseTuple(args, "ffff", &arg1, &arg2, &arg3, &arg4));
    glClearAccum(arg1, arg2, arg3, arg4);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_ClearColor(PyObject * self, PyObject * args)
{
    GLclampf arg1, arg2, arg3, arg4;
    TRY(PyArg_ParseTuple(args, "ffff", &arg1, &arg2, &arg3, &arg4));
    glClearColor(arg1, arg2, arg3, arg4);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_ClearDepth(PyObject * self, PyObject * args)
{
    GLclampd arg1;
    TRY(PyArg_ParseTuple(args, "d", &arg1));
    glClearDepth(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_ClearIndex(PyObject * self, PyObject * args)
{
    GLfloat arg1;
    TRY(PyArg_ParseTuple(args, "f", &arg1));
    glClearIndex(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_ClearStencil(PyObject * self, PyObject * args)
{
    GLint arg1;
    TRY(PyArg_ParseTuple(args, "i", &arg1));
    glClearStencil(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Color3d(PyObject * self, PyObject * args)
{
    GLdouble arg1, arg2, arg3, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "ddd", &arg1, &arg2, &arg3))
	glColor3d(arg1, arg2, arg3);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsDoubleArray(&op, &vert, &n));
	if (n < 3) {
	    PyErr_SetString(py_gl_Error, "need element with at least 3 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glColor3dv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Color3f(PyObject * self, PyObject * args)
{
    GLfloat arg1, arg2, arg3, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "fff", &arg1, &arg2, &arg3))
	glColor3f(arg1, arg2, arg3);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsFloatArray(&op, &vert, &n));
	if (n < 3) {
	    PyErr_SetString(py_gl_Error, "need element with at least 3 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glColor3fv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Color3ub(PyObject * self, PyObject * args)
{
    GLint arg1, arg2, arg3, *vert;
    GLubyte a1, a2, a3, v[3];
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "iii", &arg1, &arg2, &arg3)) {
	a1 = (GLubyte) arg1;
	a2 = (GLubyte) arg2;
	a3 = (GLubyte) arg3;
	glColor3ub(a1, a2, a3);
    } else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsIntArray(&op, &vert, &n));
	if (n < 3) {
	    PyErr_SetString(py_gl_Error, "need element with at least 3 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	v[0] = (GLubyte) vert[0];
	v[1] = (GLubyte) vert[1];
	v[2] = (GLubyte) vert[2];
	glColor3ubv(v);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Color3b(PyObject * self, PyObject * args)
{
    GLint arg1, arg2, arg3, *vert;
    GLbyte a1, a2, a3, v[3];
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "iii", &arg1, &arg2, &arg3)) {
	a1 = (GLbyte) arg1;
	a2 = (GLbyte) arg2;
	a3 = (GLbyte) arg3;
	glColor3b(a1, a2, a3);
    } else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsIntArray(&op, &vert, &n));
	if (n < 3) {
	    PyErr_SetString(py_gl_Error, "need element with at least 3 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	v[0] = (GLbyte) vert[0];
	v[1] = (GLbyte) vert[1];
	v[2] = (GLbyte) vert[2];
	glColor3bv(v);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Color3ui(PyObject * self, PyObject * args)
{
    GLint arg1, arg2, arg3, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "iii", &arg1, &arg2, &arg3)) {
	glColor3ui(arg1, arg2, arg3);
    } else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsIntArray(&op, &vert, &n));
	if (n < 3) {
	    PyErr_SetString(py_gl_Error, "need element with at least 3 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glColor3uiv((GLuint *)vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Color3i(PyObject * self, PyObject * args)
{
    GLint arg1, arg2, arg3, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "iii", &arg1, &arg2, &arg3))
	glColor3i(arg1, arg2, arg3);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsIntArray(&op, &vert, &n));
	if (n < 3) {
	    PyErr_SetString(py_gl_Error, "need element with at least 3 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glColor3iv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Color3us(PyObject * self, PyObject * args)
{
    GLshort arg1, arg2, arg3, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "hhh", &arg1, &arg2, &arg3))
	glColor3us(arg1, arg2, arg3);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsShortArray(&op, &vert, &n));
	if (n < 3) {
	    PyErr_SetString(py_gl_Error, "need element with at least 3 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glColor3usv((GLushort *)vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Color3s(PyObject * self, PyObject * args)
{
    GLshort arg1, arg2, arg3, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "hhh", &arg1, &arg2, &arg3))
	glColor3s(arg1, arg2, arg3);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsShortArray(&op, &vert, &n));
	if (n < 3) {
	    PyErr_SetString(py_gl_Error, "need element with at least 3 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glColor3sv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Color4d(PyObject * self, PyObject * args)
{
    GLdouble arg1, arg2, arg3, arg4, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "dddd", &arg1, &arg2, &arg3, &arg4))
	glColor4d(arg1, arg2, arg3, arg4);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsDoubleArray(&op, &vert, &n));
	if (n < 4) {
	    PyErr_SetString(py_gl_Error, "need element with at least 4 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glColor4dv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Color4f(PyObject * self, PyObject * args)
{
    GLfloat arg1, arg2, arg3, arg4, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "ffff", &arg1, &arg2, &arg3, &arg4))
	glColor4f(arg1, arg2, arg3, arg4);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsFloatArray(&op, &vert, &n));
	if (n < 4) {
	    PyErr_SetString(py_gl_Error, "need element with at least 4 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glColor4fv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Color4ub(PyObject * self, PyObject * args)
{
    GLint arg1, arg2, arg3, arg4, *vert;
    GLubyte a1, a2, a3, a4, v[4];
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "iiii", &arg1, &arg2, &arg3, &arg4)) {
	a1 = (GLubyte) arg1;
	a2 = (GLubyte) arg2;
	a3 = (GLubyte) arg3;
	a4 = (GLubyte) arg4;
	glColor4ub(a1, a2, a3, a4);
    } else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsIntArray(&op, &vert, &n));
	if (n < 4) {
	    PyErr_SetString(py_gl_Error, "need element with at least 4 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	v[0] = (GLubyte) arg1;
	v[1] = (GLubyte) arg2;
	v[2] = (GLubyte) arg3;
	v[3] = (GLubyte) arg4;
	glColor4ubv(v);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Color4b(PyObject * self, PyObject * args)
{
    GLint arg1, arg2, arg3, arg4, *vert;
    GLbyte a1, a2, a3, a4, v[4];
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "iiii", &arg1, &arg2, &arg3, &arg4)) {
	a1 = (GLbyte) arg1;
	a2 = (GLbyte) arg2;
	a3 = (GLbyte) arg3;
	a4 = (GLbyte) arg4;
	glColor4ub(a1, a2, a3, a4);
    } else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsIntArray(&op, &vert, &n));
	if (n < 4) {
	    PyErr_SetString(py_gl_Error, "need element with at least 4 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	v[0] = (GLbyte) arg1;
	v[1] = (GLbyte) arg2;
	v[2] = (GLbyte) arg3;
	v[3] = (GLbyte) arg4;
	glColor4bv(v);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Color4ui(PyObject * self, PyObject * args)
{
    GLint arg1, arg2, arg3, arg4, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "iiii", &arg1, &arg2, &arg3, &arg4))
	glColor4ui(arg1, arg2, arg3, arg4);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsIntArray(&op, &vert, &n));
	if (n < 4) {
	    PyErr_SetString(py_gl_Error, "need element with at least 4 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glColor4uiv((GLuint *)vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Color4i(PyObject * self, PyObject * args)
{
    GLint arg1, arg2, arg3, arg4, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "iiii", &arg1, &arg2, &arg3, &arg4))
	glColor4i(arg1, arg2, arg3, arg4);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsIntArray(&op, &vert, &n));
	if (n < 4) {
	    PyErr_SetString(py_gl_Error, "need element with at least 4 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glColor4iv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Color4us(PyObject * self, PyObject * args)
{
    GLshort arg1, arg2, arg3, arg4, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "hhhh", &arg1, &arg2, &arg3, &arg4))
	glColor4us(arg1, arg2, arg3, arg4);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsShortArray(&op, &vert, &n));
	if (n < 4) {
	    PyErr_SetString(py_gl_Error, "need element with at least 4 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glColor4usv((GLushort *)vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Color4s(PyObject * self, PyObject * args)
{
    GLshort arg1, arg2, arg3, arg4, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "hhhh", &arg1, &arg2, &arg3, &arg4))
	glColor4s(arg1, arg2, arg3, arg4);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsShortArray(&op, &vert, &n));
	if (n < 4) {
	    PyErr_SetString(py_gl_Error, "need element with at least 4 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glColor4sv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_ColorMask(PyObject * self, PyObject * args)
{
    GLboolean arg1, arg2, arg3, arg4;
    TRY(PyArg_ParseTuple(args, "bbbb", &arg1, &arg2, &arg3, &arg4));
    glColorMask(arg1, arg2, arg3, arg4);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_ColorMaterial(PyObject * self, PyObject * args)
{
    GLenum arg1, arg2;
    TRY(PyArg_ParseTuple(args, "ii", &arg1, &arg2));
    glColorMaterial(arg1, arg2);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_CopyPixels(PyObject * self, PyObject * args)
{
    GLint arg1, arg2;
    GLsizei arg3, arg4;
    GLenum arg5;
    TRY(PyArg_ParseTuple(args, "iiiii", &arg1, &arg2, &arg3, &arg4, &arg5));
    glCopyPixels(arg1, arg2, arg3, arg4, arg5);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_CullFace(PyObject * self, PyObject * args)
{
    GLenum arg1;
    TRY(PyArg_ParseTuple(args, "i", &arg1));
    glCullFace(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_DeleteLists(PyObject * self, PyObject * args)
{
    GLuint arg1;
    GLsizei arg2;
    TRY(PyArg_ParseTuple(args, "ii", &arg1, &arg2));
    glDeleteLists(arg1, arg2);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_DepthFunc(PyObject * self, PyObject * args)
{
    GLenum arg1;
    TRY(PyArg_ParseTuple(args, "i", &arg1));
    glDepthFunc(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_DepthMask(PyObject * self, PyObject * args)
{
    GLboolean arg1;
    TRY(PyArg_ParseTuple(args, "b", &arg1));
    glDepthMask(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_DepthRange(PyObject * self, PyObject * args)
{
    GLclampd arg1, arg2;
    TRY(PyArg_ParseTuple(args, "dd", &arg1, &arg2));
    glDepthRange(arg1, arg2);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Disable(PyObject * self, PyObject * args)
{
    GLenum arg1;
    TRY(PyArg_ParseTuple(args, "i", &arg1));
    glDisable(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_DrawBuffer(PyObject * self, PyObject * args)
{
    GLenum arg1;
    TRY(PyArg_ParseTuple(args, "i", &arg1));
    glDrawBuffer(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_EdgeFlag(PyObject * self, PyObject * args)
{
    GLboolean arg1;
    TRY(PyArg_ParseTuple(args, "b", &arg1));
    glEdgeFlag(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Enable(PyObject * self, PyObject * args)
{
    GLenum arg1;
    TRY(PyArg_ParseTuple(args, "i", &arg1));
    glEnable(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_End(PyObject * self, PyObject * args)
{
    TRY(PyArg_ParseTuple(args, ""));
    glEnd();
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_EndList(PyObject * self, PyObject * args)
{
    TRY(PyArg_ParseTuple(args, ""));
    glEndList();
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_EvalCoord1d(PyObject * self, PyObject * args)
{
    GLdouble arg1;
    TRY(PyArg_ParseTuple(args, "d", &arg1));
    glEvalCoord1d(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_EvalCoord1f(PyObject * self, PyObject * args)
{
    GLfloat arg1;
    TRY(PyArg_ParseTuple(args, "f", &arg1));
    glEvalCoord1f(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_EvalCoord2d(PyObject * self, PyObject * args)
{
    GLdouble arg1, arg2, *coords;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "dd", &arg1, &arg2))
	glEvalCoord2d(arg1, arg2);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsDoubleArray(&op, &coords, &n));
	if (n < 2) {
	    PyErr_SetString(py_gl_Error, "need element with at least 2 items");
	    PyArray_ClearMemory(op, coords);
	    return NULL;
	}
	glEvalCoord2dv(coords);
	PyArray_ClearMemory(op, coords);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_EvalCoord2f(PyObject * self, PyObject * args)
{
    GLfloat arg1, arg2, *coords;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "ff", &arg1, &arg2))
	glEvalCoord2f(arg1, arg2);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsFloatArray(&op, &coords, &n));
	if (n < 2) {
	    PyErr_SetString(py_gl_Error, "need element with at least 2 items");
	    PyArray_ClearMemory(op, coords);
	    return NULL;
	}
	glEvalCoord2fv(coords);
	PyArray_ClearMemory(op, coords);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_EvalMesh1(PyObject * self, PyObject * args)
{
    GLenum arg1;
    GLint arg2, arg3;
    TRY(PyArg_ParseTuple(args, "iii", &arg1, &arg2, &arg3));
    glEvalMesh1(arg1, arg2, arg3);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_EvalMesh2(PyObject * self, PyObject * args)
{
    GLenum arg1;
    GLint arg2, arg3, arg4, arg5;
    TRY(PyArg_ParseTuple(args, "iiiii", &arg1, &arg2, &arg3, &arg4, &arg5));
    glEvalMesh2(arg1, arg2, arg3, arg4, arg5);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_EvalPoint1(PyObject * self, PyObject * args)
{
    GLint arg1;
    TRY(PyArg_ParseTuple(args, "i", &arg1));
    glEvalPoint1(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_EvalPoint2(PyObject * self, PyObject * args)
{
    GLint arg1, arg2;
    TRY(PyArg_ParseTuple(args, "ii", &arg1, &arg2));
    glEvalPoint2(arg1, arg2);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Finish(PyObject * self, PyObject * args)
{
    TRY(PyArg_ParseTuple(args, ""));
    glFinish();
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Flush(PyObject * self, PyObject * args)
{
    TRY(PyArg_ParseTuple(args, ""));
    glFlush();
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Fogf(PyObject * self, PyObject * args)
{
    GLenum arg1;
    GLfloat arg2, *rest;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "if", &arg1, &arg2))
	glFogf(arg1, arg2);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "iO", &arg1, &op));
	TRY(PyArray_AsFloatArray(&op, &rest, &n));
	if ((n == 1) || (n > 3)) {
	    glFogfv(arg1, rest);
	    PyArray_ClearMemory(op, rest);
	} else {
	    PyErr_SetString(py_gl_Error,
			    "second argument needs 1 or at least 4 items!");
	    PyArray_ClearMemory(op, rest);
	    return NULL;
	}
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Fogi(PyObject * self, PyObject * args)
{
    GLenum arg1;
    GLint arg2, *rest;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "ii", &arg1, &arg2))
	glFogi(arg1, arg2);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "iO", &arg1, &op));
	TRY(PyArray_AsIntArray(&op, &rest, &n));
	if ((n == 1) || (n > 3)) {
	    glFogiv(arg1, rest);
	    PyArray_ClearMemory(op, rest);
	} else {
	    PyErr_SetString(py_gl_Error,
			    "second argument needs 1 or at least 4 items!");
	    PyArray_ClearMemory(op, rest);
	    return NULL;
	}
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_FrontFace(PyObject * self, PyObject * args)
{
    GLenum arg1;
    TRY(PyArg_ParseTuple(args, "i", &arg1));
    glFrontFace(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Frustum(PyObject * self, PyObject * args)
{
    GLdouble arg1, arg2, arg3, arg4, arg5, arg6;
    TRY(PyArg_ParseTuple
	(args, "dddddd", &arg1, &arg2, &arg3, &arg4, &arg5, &arg6));
    glFrustum(arg1, arg2, arg3, arg4, arg5, arg6);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_GetString(PyObject * self, PyObject * args)
{
    /* WARNING: only returns values after an OpenGL context has been created */
    GLenum arg1;
    char *p;

    TRY(PyArg_ParseTuple(args, "i", &arg1));
    p = (char *) glGetString(arg1);
    return Py_BuildValue("s", p);
}

static PyObject *py_gl_Hint(PyObject * self, PyObject * args)
{
    GLenum arg1, arg2;
    TRY(PyArg_ParseTuple(args, "ii", &arg1, &arg2));
    glHint(arg1, arg2);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_IndexMask(PyObject * self, PyObject * args)
{
    GLuint arg1;
    TRY(PyArg_ParseTuple(args, "i", &arg1));
    glIndexMask(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Indexd(PyObject * self, PyObject * args)
{
    GLdouble arg1;
    TRY(PyArg_ParseTuple(args, "d", &arg1));
    glIndexd(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Indexf(PyObject * self, PyObject * args)
{
    GLfloat arg1;
    TRY(PyArg_ParseTuple(args, "f", &arg1));
    glIndexf(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Indexi(PyObject * self, PyObject * args)
{
    GLint arg1;
    TRY(PyArg_ParseTuple(args, "i", &arg1));
    glIndexi(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Indexs(PyObject * self, PyObject * args)
{
    GLshort arg1;
    TRY(PyArg_ParseTuple(args, "h", &arg1));
    glIndexs(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_InitNames(PyObject * self, PyObject * args)
{
    TRY(PyArg_ParseTuple(args, ""));
    glInitNames();
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_LightModelf(PyObject * self, PyObject * args)
{
    GLenum arg1;
    GLfloat arg2, *rest;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "if", &arg1, &arg2))
	glLightModelf(arg1, arg2);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "iO", &arg1, &op));
	TRY(PyArray_AsFloatArray(&op, &rest, &n));
	if ((n == 1) || (n > 3)) {
	    glLightModelfv(arg1, rest);
	    PyArray_ClearMemory(op, rest);
	} else {
	    PyErr_SetString(py_gl_Error,
			    "second argument needs 1 or at least 4 items!");
	    PyArray_ClearMemory(op, rest);
	    return NULL;
	}
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_LightModeli(PyObject * self, PyObject * args)
{
    GLenum arg1;
    GLint arg2, *rest;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "ii", &arg1, &arg2))
	glLightModeli(arg1, arg2);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "iO", &arg1, &op));
	TRY(PyArray_AsIntArray(&op, &rest, &n));
	if ((n == 1) || (n > 3)) {
	    glLightModeliv(arg1, rest);
	    PyArray_ClearMemory(op, rest);
	} else {
	    PyErr_SetString(py_gl_Error,
			    "second argument needs 1 or at least 4 items!");
	    PyArray_ClearMemory(op, rest);
	    return NULL;
	}
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Lightf(PyObject * self, PyObject * args)
{
    GLenum arg1, arg2;
    GLfloat arg3, *rest;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "iif", &arg1, &arg2, &arg3))
	glLightf(arg1, arg2, arg3);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "iiO", &arg1, &arg2, &op));
	TRY(PyArray_AsFloatArray(&op, &rest, &n));
	if ((n == 1) || (n > 2)) {
	    glLightfv(arg1, arg2, rest);
	    PyArray_ClearMemory(op, rest);
	} else {
	    PyErr_SetString(py_gl_Error,
			    "second argument needs 1 or at least 3 items!");
	    PyArray_ClearMemory(op, rest);
	    return NULL;
	}
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Lighti(PyObject * self, PyObject * args)
{
    GLenum arg1, arg2;
    GLint arg3, *rest;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "iii", &arg1, &arg2, &arg3))
	glLighti(arg1, arg2, arg3);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "iiO", &arg1, &arg2, &op));
	TRY(PyArray_AsIntArray(&op, &rest, &n));
	if ((n == 1) || (n > 2)) {
	    glLightiv(arg1, arg2, rest);
	    PyArray_ClearMemory(op, rest);
	} else {
	    PyErr_SetString(py_gl_Error,
			    "second argument needs 1 or at least 3 items!");
	    PyArray_ClearMemory(op, rest);
	    return NULL;
	}
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_LineStipple(PyObject * self, PyObject * args)
{
    GLint arg1;
    GLushort arg2;
    TRY(PyArg_ParseTuple(args, "ih", &arg1, &arg2));
    glLineStipple(arg1, arg2);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_LineWidth(PyObject * self, PyObject * args)
{
    GLfloat arg1;
    TRY(PyArg_ParseTuple(args, "f", &arg1));
    glLineWidth(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_ListBase(PyObject * self, PyObject * args)
{
    GLuint arg1;
    TRY(PyArg_ParseTuple(args, "i", &arg1));
    glListBase(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_LoadIdentity(PyObject * self, PyObject * args)
{
    TRY(PyArg_ParseTuple(args, ""));
    glLoadIdentity();
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_LoadName(PyObject * self, PyObject * args)
{
    GLuint arg1;
    TRY(PyArg_ParseTuple(args, "i", &arg1));
    glLoadName(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_LogicOp(PyObject * self, PyObject * args)
{
    GLenum arg1;
    TRY(PyArg_ParseTuple(args, "i", &arg1));
    glLogicOp(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *py_gl_Map1d(PyObject * self, PyObject * args)
{

    GLenum arg1;
    GLdouble arg2;
    GLdouble arg3;
    GLint arg4;
    GLint arg5;
    PyObject *arg6;
    GLdouble *points;
    int n;

    TRY(PyArg_ParseTuple
	(args, "iddiiO", &arg1, &arg2, &arg3, &arg4, &arg5, &arg6));
    TRY(PyArray_AsDoubleArray(&arg6, &points, &n));
    glMap1d(arg1, arg2, arg3, arg4, arg5, points);
    Py_DECREF(arg6);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Map1f(PyObject * self, PyObject * args)
{
    GLenum arg1;
    GLfloat arg2;
    GLfloat arg3;
    GLint arg4;
    GLint arg5;
    PyObject *arg6;
    GLfloat *points;
    int n;

    TRY(PyArg_ParseTuple
	(args, "iffiiO", &arg1, &arg2, &arg3, &arg4, &arg5, &arg6));
    TRY(PyArray_AsFloatArray(&arg6, &points, &n));
    glMap1f(arg1, arg2, arg3, arg4, arg5, points);
    Py_DECREF(arg6);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Map2d(PyObject * self, PyObject * args)
{
    GLenum arg1;
    GLdouble arg2;
    GLdouble arg3;
    GLint arg4;
    GLint arg5;
    GLdouble arg6;
    GLdouble arg7;
    GLint arg8;
    GLint arg9;
    PyObject *arg10;
    GLdouble *points;
    int n;

    TRY(PyArg_ParseTuple(args, "iddiiddiiO",
			 &arg1, &arg2, &arg3, &arg4, &arg5,
			 &arg6, &arg7, &arg8, &arg9, &arg10));

    TRY(PyArray_AsDoubleArray(&arg10, &points, &n));
    glMap2d(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, points);

    Py_DECREF(arg10);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Map2f(PyObject * self, PyObject * args)
{
    GLenum arg1;
    GLfloat arg2;
    GLfloat arg3;
    GLint arg4;
    GLint arg5;
    GLfloat arg6;
    GLfloat arg7;
    GLint arg8;
    GLint arg9;
    PyObject *arg10;
    GLfloat *points;
    int n;

    TRY(PyArg_ParseTuple(args, "iffiiffiiO",
			 &arg1, &arg2, &arg3, &arg4, &arg5,
			 &arg6, &arg7, &arg8, &arg9, &arg10));

    TRY(PyArray_AsFloatArray(&arg10, &points, &n));
    glMap2f(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, points);

    Py_DECREF(arg10);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_MapGrid1d(PyObject * self, PyObject * args)
{
    GLint arg1;
    GLdouble arg2, arg3;
    TRY(PyArg_ParseTuple(args, "idd", &arg1, &arg2, &arg3));
    glMapGrid1d(arg1, arg2, arg3);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_MapGrid1f(PyObject * self, PyObject * args)
{
    GLint arg1;
    GLfloat arg2, arg3;
    TRY(PyArg_ParseTuple(args, "iff", &arg1, &arg2, &arg3));
    glMapGrid1f(arg1, arg2, arg3);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_MapGrid2d(PyObject * self, PyObject * args)
{
    GLint arg1, arg4;
    GLdouble arg2, arg3, arg5, arg6;
    TRY(PyArg_ParseTuple
	(args, "iddidd", &arg1, &arg2, &arg3, &arg4, &arg5, &arg6));
    glMapGrid2d(arg1, arg2, arg3, arg4, arg5, arg6);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_MapGrid2f(PyObject * self, PyObject * args)
{
    GLint arg1, arg4;
    GLfloat arg2, arg3, arg5, arg6;
    TRY(PyArg_ParseTuple
	(args, "iffiff", &arg1, &arg2, &arg3, &arg4, &arg5, &arg6));
    glMapGrid2f(arg1, arg2, arg3, arg4, arg5, arg6);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Materialf(PyObject * self, PyObject * args)
{
    GLenum arg1, arg2;
    GLfloat arg3, *rest;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "iif", &arg1, &arg2, &arg3))
	glMaterialf(arg1, arg2, arg3);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "iiO", &arg1, &arg2, &op));
	TRY(PyArray_AsFloatArray(&op, &rest, &n));
	if ((n == 1) || (n > 3)) {
	    glMaterialfv(arg1, arg2, rest);
	    PyArray_ClearMemory(op, rest);
	} else {
	    PyErr_SetString(py_gl_Error,
			    "second argument needs 1 or at least 4 items!");
	    PyArray_ClearMemory(op, rest);
	    return NULL;
	}
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Materiali(PyObject * self, PyObject * args)
{
    GLenum arg1, arg2;
    GLint arg3, *rest;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "iii", &arg1, &arg2, &arg3))
	glMateriali(arg1, arg2, arg3);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "iiO", &arg1, &arg2, &op));
	TRY(PyArray_AsIntArray(&op, &rest, &n));
	if ((n == 1) || (n > 3)) {
	    glMaterialiv(arg1, arg2, rest);
	    PyArray_ClearMemory(op, rest);
	} else {
	    PyErr_SetString(py_gl_Error,
			    "second argument needs 1 or at least 4 items!");
	    PyArray_ClearMemory(op, rest);
	    return NULL;
	}
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_MatrixMode(PyObject * self, PyObject * args)
{
    GLenum arg1;
    TRY(PyArg_ParseTuple(args, "i", &arg1));
    glMatrixMode(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_MultMatrixd(PyObject * self, PyObject * args)
{
    GLdouble *mat;
    PyObject *op;
    int n;
    TRY(PyArg_ParseTuple(args, "O", &op));
    TRY(PyArray_AsDoubleArray(&op, &mat, &n));
    if (n < 16) {
	PyErr_SetString(py_gl_Error, "need element with at least 16 items");
	PyArray_ClearMemory(op, mat);
	return NULL;
    }
    glMultMatrixd(mat);
    PyArray_ClearMemory(op, mat);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_MultMatrixf(PyObject * self, PyObject * args)
{
    GLfloat *mat;
    PyObject *op;
    int n;
    TRY(PyArg_ParseTuple(args, "O", &op));
    TRY(PyArray_AsFloatArray(&op, &mat, &n));
    if (n < 16) {
	PyErr_SetString(py_gl_Error, "need element with at least 16 items");
	PyArray_ClearMemory(op, mat);
	return NULL;
    }
    glMultMatrixf(mat);
    PyArray_ClearMemory(op, mat);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_NewList(PyObject * self, PyObject * args)
{
    GLuint arg1;
    GLenum arg2;
    TRY(PyArg_ParseTuple(args, "ii", &arg1, &arg2));
    glNewList(arg1, arg2);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Normal3d(PyObject * self, PyObject * args)
{
    GLdouble arg1, arg2, arg3, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "ddd", &arg1, &arg2, &arg3))
	glNormal3d(arg1, arg2, arg3);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsDoubleArray(&op, &vert, &n));
	if (n < 3) {
	    PyErr_SetString(py_gl_Error, "need element with at least 3 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glNormal3dv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Normal3f(PyObject * self, PyObject * args)
{
    GLfloat arg1, arg2, arg3, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "fff", &arg1, &arg2, &arg3))
	glNormal3f(arg1, arg2, arg3);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsFloatArray(&op, &vert, &n));
	if (n < 3) {
	    PyErr_SetString(py_gl_Error, "need element with at least 3 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glNormal3fv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Normal3i(PyObject * self, PyObject * args)
{
    GLint arg1, arg2, arg3, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "iii", &arg1, &arg2, &arg3))
	glNormal3i(arg1, arg2, arg3);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsIntArray(&op, &vert, &n));
	if (n < 3) {
	    PyErr_SetString(py_gl_Error, "need element with at least 3 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glNormal3iv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Normal3s(PyObject * self, PyObject * args)
{
    GLshort arg1, arg2, arg3, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "hhh", &arg1, &arg2, &arg3))
	glNormal3s(arg1, arg2, arg3);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsShortArray(&op, &vert, &n));
	if (n < 3) {
	    PyErr_SetString(py_gl_Error, "need element with at least 3 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glNormal3sv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Ortho(PyObject * self, PyObject * args)
{
    GLdouble arg1, arg2, arg3, arg4, arg5, arg6;
    TRY(PyArg_ParseTuple
	(args, "dddddd", &arg1, &arg2, &arg3, &arg4, &arg5, &arg6));
    glOrtho(arg1, arg2, arg3, arg4, arg5, arg6);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_PassThrough(PyObject * self, PyObject * args)
{
    GLfloat arg1;
    TRY(PyArg_ParseTuple(args, "f", &arg1));
    glPassThrough(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_PixelStoref(PyObject * self, PyObject * args)
{
    GLenum arg1;
    GLfloat arg2;
    TRY(PyArg_ParseTuple(args, "if", &arg1, &arg2));
    glPixelStoref(arg1, arg2);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_PixelStorei(PyObject * self, PyObject * args)
{
    GLenum arg1;
    GLint arg2;
    TRY(PyArg_ParseTuple(args, "ii", &arg1, &arg2));
    glPixelStorei(arg1, arg2);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_PixelTransferf(PyObject * self, PyObject * args)
{
    GLenum arg1;
    GLfloat arg2;
    TRY(PyArg_ParseTuple(args, "if", &arg1, &arg2));
    glPixelTransferf(arg1, arg2);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_PixelTransferi(PyObject * self, PyObject * args)
{
    GLenum arg1;
    GLint arg2;
    TRY(PyArg_ParseTuple(args, "ii", &arg1, &arg2));
    glPixelTransferi(arg1, arg2);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_PixelZoom(PyObject * self, PyObject * args)
{
    GLfloat arg1, arg2;
    TRY(PyArg_ParseTuple(args, "ff", &arg1, &arg2));
    glPixelZoom(arg1, arg2);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_PointSize(PyObject * self, PyObject * args)
{
    GLfloat arg1;
    TRY(PyArg_ParseTuple(args, "f", &arg1));
    glPointSize(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

#ifdef NUMERIC
static PyObject *py_gl_PolygonStipple(PyObject * self, PyObject * args)
{
    PyObject *op;
    GLubyte *mask;
    int n;

    TRY(PyArg_ParseTuple(args, "O", &op));
    TRY(PyArray_AsUByteArray(&op, &mask, &n));
    glPolygonStipple(mask);
    Py_INCREF(Py_None);
    return Py_None;
}

#endif

static PyObject *py_gl_PolygonMode(PyObject * self, PyObject * args)
{
    GLenum arg1, arg2;
    TRY(PyArg_ParseTuple(args, "ii", &arg1, &arg2));
    glPolygonMode(arg1, arg2);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_PopAttrib(PyObject * self, PyObject * args)
{
    TRY(PyArg_ParseTuple(args, ""));
    glPopAttrib();
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_PopMatrix(PyObject * self, PyObject * args)
{
    TRY(PyArg_ParseTuple(args, ""));
    glPopMatrix();
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_PopName(PyObject * self, PyObject * args)
{
    TRY(PyArg_ParseTuple(args, ""));
    glPopName();
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_PushAttrib(PyObject * self, PyObject * args)
{
    GLbitfield arg1;
    TRY(PyArg_ParseTuple(args, "i", &arg1));
    glPushAttrib(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_PushMatrix(PyObject * self, PyObject * args)
{
    TRY(PyArg_ParseTuple(args, ""));
    glPushMatrix();
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_PushName(PyObject * self, PyObject * args)
{
    GLuint arg1;
    TRY(PyArg_ParseTuple(args, "i", &arg1));
    glPushName(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_RasterPos2d(PyObject * self, PyObject * args)
{
    GLdouble arg1, arg2, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "dd", &arg1, &arg2))
	glRasterPos2d(arg1, arg2);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsDoubleArray(&op, &vert, &n));
	if (n < 2) {
	    PyErr_SetString(py_gl_Error, "need element with at least 2 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glRasterPos2dv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_RasterPos2f(PyObject * self, PyObject * args)
{
    GLfloat arg1, arg2, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "ff", &arg1, &arg2))
	glRasterPos2f(arg1, arg2);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsFloatArray(&op, &vert, &n));
	if (n < 2) {
	    PyErr_SetString(py_gl_Error, "need element with at least 2 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glRasterPos2fv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_RasterPos2i(PyObject * self, PyObject * args)
{
    GLint arg1, arg2, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "ii", &arg1, &arg2))
	glRasterPos2i(arg1, arg2);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsIntArray(&op, &vert, &n));
	if (n < 2) {
	    PyErr_SetString(py_gl_Error, "need element with at least 2 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glRasterPos2iv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_RasterPos2s(PyObject * self, PyObject * args)
{
    GLshort arg1, arg2, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "hh", &arg1, &arg2))
	glRasterPos2s(arg1, arg2);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsShortArray(&op, &vert, &n));
	if (n < 2) {
	    PyErr_SetString(py_gl_Error, "need element with at least 2 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glRasterPos2sv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_RasterPos3d(PyObject * self, PyObject * args)
{
    GLdouble arg1, arg2, arg3, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "ddd", &arg1, &arg2, &arg3))
	glRasterPos3d(arg1, arg2, arg3);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsDoubleArray(&op, &vert, &n));
	if (n < 3) {
	    PyErr_SetString(py_gl_Error, "need element with at least 3 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glRasterPos3dv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_RasterPos3f(PyObject * self, PyObject * args)
{
    GLfloat arg1, arg2, arg3, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "fff", &arg1, &arg2, &arg3))
	glRasterPos3f(arg1, arg2, arg3);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsFloatArray(&op, &vert, &n));
	if (n < 3) {
	    PyErr_SetString(py_gl_Error, "need element with at least 3 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glRasterPos3fv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_RasterPos3i(PyObject * self, PyObject * args)
{
    GLint arg1, arg2, arg3, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "iii", &arg1, &arg2, &arg3))
	glRasterPos3i(arg1, arg2, arg3);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsIntArray(&op, &vert, &n));
	if (n < 3) {
	    PyErr_SetString(py_gl_Error, "need element with at least 3 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glRasterPos3iv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_RasterPos3s(PyObject * self, PyObject * args)
{
    GLshort arg1, arg2, arg3, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "hhh", &arg1, &arg2, &arg3))
	glRasterPos3s(arg1, arg2, arg3);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsShortArray(&op, &vert, &n));
	if (n < 3) {
	    PyErr_SetString(py_gl_Error, "need element with at least 3 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glRasterPos3sv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_RasterPos4d(PyObject * self, PyObject * args)
{
    GLdouble arg1, arg2, arg3, arg4, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "dddd", &arg1, &arg2, &arg3, &arg4))
	glRasterPos4d(arg1, arg2, arg3, arg4);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsDoubleArray(&op, &vert, &n));
	if (n < 4) {
	    PyErr_SetString(py_gl_Error, "need element with at least 4 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glRasterPos4dv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_RasterPos4f(PyObject * self, PyObject * args)
{
    GLfloat arg1, arg2, arg3, arg4, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "ffff", &arg1, &arg2, &arg3, &arg4))
	glRasterPos4f(arg1, arg2, arg3, arg4);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsFloatArray(&op, &vert, &n));
	if (n < 4) {
	    PyErr_SetString(py_gl_Error, "need element with at least 4 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glRasterPos4fv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_RasterPos4i(PyObject * self, PyObject * args)
{
    GLint arg1, arg2, arg3, arg4, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "iiii", &arg1, &arg2, &arg3, &arg4))
	glRasterPos4i(arg1, arg2, arg3, arg4);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsIntArray(&op, &vert, &n));
	if (n < 4) {
	    PyErr_SetString(py_gl_Error, "need element with at least 4 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glRasterPos4iv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_RasterPos4s(PyObject * self, PyObject * args)
{
    GLshort arg1, arg2, arg3, arg4, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "hhhh", &arg1, &arg2, &arg3, &arg4))
	glRasterPos4s(arg1, arg2, arg3, arg4);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsShortArray(&op, &vert, &n));
	if (n < 4) {
	    PyErr_SetString(py_gl_Error, "need element with at least 4 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glRasterPos4sv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_ReadBuffer(PyObject * self, PyObject * args)
{
    GLenum arg1;
    TRY(PyArg_ParseTuple(args, "i", &arg1));
    glReadBuffer(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Rectd(PyObject * self, PyObject * args)
{
    GLdouble arg1, arg2, arg3, arg4, *vert1, *vert2;
    PyObject *op1, *op2;
    int n1, n2;
    if (PyArg_ParseTuple(args, "dddd", &arg1, &arg2, &arg3, &arg4))
	glRectd(arg1, arg2, arg3, arg4);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "OO", &op1, &op2));
	TRY(PyArray_AsDoubleArray(&op1, &vert1, &n1));
	TRY(PyArray_AsDoubleArray(&op2, &vert2, &n2));
	if ((n1 < 2) || (n2 < 2)) {
	    PyErr_SetString(py_gl_Error, "need elements with at least 2 items");
	    Py_DECREF(op1);
	    Py_DECREF(op2);
	    return NULL;
	}
	glRectdv(vert1, vert2);
	Py_DECREF(op1);
	Py_DECREF(op2);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Rectf(PyObject * self, PyObject * args)
{
    GLfloat arg1, arg2, arg3, arg4, *vert1, *vert2;
    PyObject *op1, *op2;
    int n1, n2;
    if (PyArg_ParseTuple(args, "ffff", &arg1, &arg2, &arg3, &arg4))
	glRectf(arg1, arg2, arg3, arg4);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "OO", &op1, &op2));
	TRY(PyArray_AsFloatArray(&op1, &vert1, &n1));
	TRY(PyArray_AsFloatArray(&op2, &vert2, &n2));
	if ((n1 < 2) || (n2 < 2)) {
	    PyErr_SetString(py_gl_Error, "need elements with at least 2 items");
	    Py_DECREF(op1);
	    Py_DECREF(op2);
	    return NULL;
	}
	glRectfv(vert1, vert2);
	Py_DECREF(op1);
	Py_DECREF(op2);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Recti(PyObject * self, PyObject * args)
{
    GLint arg1, arg2, arg3, arg4, *vert1, *vert2;
    PyObject *op1, *op2;
    int n1, n2;
    if (PyArg_ParseTuple(args, "iiii", &arg1, &arg2, &arg3, &arg4))
	glRecti(arg1, arg2, arg3, arg4);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "OO", &op1, &op2));
	TRY(PyArray_AsIntArray(&op1, &vert1, &n1));
	TRY(PyArray_AsIntArray(&op2, &vert2, &n2));
	if ((n1 < 2) || (n2 < 2)) {
	    PyErr_SetString(py_gl_Error, "need elements with at least 2 items");
	    Py_DECREF(op1);
	    Py_DECREF(op2);
	    return NULL;
	}
	glRectiv(vert1, vert2);
	Py_DECREF(op1);
	Py_DECREF(op2);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Rects(PyObject * self, PyObject * args)
{
    GLshort arg1, arg2, arg3, arg4, *vert1, *vert2;
    PyObject *op1, *op2;
    int n1, n2;
    if (PyArg_ParseTuple(args, "hhhh", &arg1, &arg2, &arg3, &arg4))
	glRects(arg1, arg2, arg3, arg4);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "OO", &op1, &op2));
	TRY(PyArray_AsShortArray(&op1, &vert1, &n1));
	TRY(PyArray_AsShortArray(&op2, &vert2, &n2));
	if ((n1 < 2) || (n2 < 2)) {
	    PyErr_SetString(py_gl_Error, "need elements with at least 2 items");
	    Py_DECREF(op1);
	    Py_DECREF(op2);
	    return NULL;
	}
	glRectsv(vert1, vert2);
	Py_DECREF(op1);
	Py_DECREF(op2);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Rotated(PyObject * self, PyObject * args)
{
    GLdouble arg1, arg2, arg3, arg4;
    TRY(PyArg_ParseTuple(args, "dddd", &arg1, &arg2, &arg3, &arg4));
    glRotated(arg1, arg2, arg3, arg4);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Rotatef(PyObject * self, PyObject * args)
{
    GLfloat arg1, arg2, arg3, arg4;
    TRY(PyArg_ParseTuple(args, "ffff", &arg1, &arg2, &arg3, &arg4));
    glRotatef(arg1, arg2, arg3, arg4);
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *py_gl_Scaled(PyObject * self, PyObject * args)
{
    GLdouble arg1, arg2, arg3;
    TRY(PyArg_ParseTuple(args, "ddd", &arg1, &arg2, &arg3));
    glScaled(arg1, arg2, arg3);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Scalef(PyObject * self, PyObject * args)
{
    GLfloat arg1, arg2, arg3;
    TRY(PyArg_ParseTuple(args, "fff", &arg1, &arg2, &arg3));
    glScalef(arg1, arg2, arg3);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Scissor(PyObject * self, PyObject * args)
{
    GLint arg1, arg2;
    GLsizei arg3, arg4;
    TRY(PyArg_ParseTuple(args, "iiii", &arg1, &arg2, &arg3, &arg4));
    glScissor(arg1, arg2, arg3, arg4);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_ShadeModel(PyObject * self, PyObject * args)
{
    GLenum arg1;
    TRY(PyArg_ParseTuple(args, "i", &arg1));
    glShadeModel(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_StencilFunc(PyObject * self, PyObject * args)
{
    GLenum arg1;
    GLint arg2;
    GLuint arg3;
    TRY(PyArg_ParseTuple(args, "iii", &arg1, &arg2, &arg3));
    glStencilFunc(arg1, arg2, arg3);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_StencilMask(PyObject * self, PyObject * args)
{
    GLuint arg1;
    TRY(PyArg_ParseTuple(args, "i", &arg1));
    glStencilMask(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_StencilOp(PyObject * self, PyObject * args)
{
    GLenum arg1, arg2, arg3;
    TRY(PyArg_ParseTuple(args, "iii", &arg1, &arg2, &arg3));
    glStencilOp(arg1, arg2, arg3);
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *py_gl_TexCoord1d(PyObject * self, PyObject * args)
{
    GLdouble arg1;
    TRY(PyArg_ParseTuple(args, "d", &arg1));
    glTexCoord1d(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_TexCoord1f(PyObject * self, PyObject * args)
{
    GLfloat arg1;
    TRY(PyArg_ParseTuple(args, "f", &arg1));
    glTexCoord1f(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_TexCoord1i(PyObject * self, PyObject * args)
{
    GLint arg1;
    TRY(PyArg_ParseTuple(args, "i", &arg1));
    glTexCoord1i(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_TexCoord1s(PyObject * self, PyObject * args)
{
    GLshort arg1;
    TRY(PyArg_ParseTuple(args, "h", &arg1));
    glTexCoord1s(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_TexCoord2d(PyObject * self, PyObject * args)
{
    GLdouble arg1, arg2, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "dd", &arg1, &arg2))
	glTexCoord2d(arg1, arg2);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsDoubleArray(&op, &vert, &n));
	if (n < 2) {
	    PyErr_SetString(py_gl_Error, "need element with at least 2 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glTexCoord2dv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_TexCoord2f(PyObject * self, PyObject * args)
{
    GLfloat arg1, arg2, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "ff", &arg1, &arg2))
	glTexCoord2f(arg1, arg2);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsFloatArray(&op, &vert, &n));
	if (n < 2) {
	    PyErr_SetString(py_gl_Error, "need element with at least 2 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glTexCoord2fv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_TexCoord2i(PyObject * self, PyObject * args)
{
    GLint arg1, arg2, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "ii", &arg1, &arg2))
	glTexCoord2i(arg1, arg2);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsIntArray(&op, &vert, &n));
	if (n < 2) {
	    PyErr_SetString(py_gl_Error, "need element with at least 2 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glTexCoord2iv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_TexCoord2s(PyObject * self, PyObject * args)
{
    GLshort arg1, arg2, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "hh", &arg1, &arg2))
	glTexCoord2s(arg1, arg2);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsShortArray(&op, &vert, &n));
	if (n < 2) {
	    PyErr_SetString(py_gl_Error, "need element with at least 2 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glTexCoord2sv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_TexCoord3d(PyObject * self, PyObject * args)
{
    GLdouble arg1, arg2, arg3, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "ddd", &arg1, &arg2, &arg3))
	glTexCoord3d(arg1, arg2, arg3);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsDoubleArray(&op, &vert, &n));
	if (n < 3) {
	    PyErr_SetString(py_gl_Error, "need element with at least 3 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glTexCoord3dv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_TexCoord3f(PyObject * self, PyObject * args)
{
    GLfloat arg1, arg2, arg3, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "fff", &arg1, &arg2, &arg3))
	glTexCoord3f(arg1, arg2, arg3);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsFloatArray(&op, &vert, &n));
	if (n < 3) {
	    PyErr_SetString(py_gl_Error, "need element with at least 3 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glTexCoord3fv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_TexCoord3i(PyObject * self, PyObject * args)
{
    GLint arg1, arg2, arg3, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "iii", &arg1, &arg2, &arg3))
	glTexCoord3i(arg1, arg2, arg3);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsIntArray(&op, &vert, &n));
	if (n < 3) {
	    PyErr_SetString(py_gl_Error, "need element with at least 3 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glTexCoord3iv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_TexCoord3s(PyObject * self, PyObject * args)
{
    GLshort arg1, arg2, arg3, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "hhh", &arg1, &arg2, &arg3))
	glTexCoord3s(arg1, arg2, arg3);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsShortArray(&op, &vert, &n));
	if (n < 3) {
	    PyErr_SetString(py_gl_Error, "need element with at least 3 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glTexCoord3sv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_TexCoord4d(PyObject * self, PyObject * args)
{
    GLdouble arg1, arg2, arg3, arg4, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "dddd", &arg1, &arg2, &arg3, &arg4))
	glTexCoord4d(arg1, arg2, arg3, arg4);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsDoubleArray(&op, &vert, &n));
	if (n < 4) {
	    PyErr_SetString(py_gl_Error, "need element with at least 4 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glTexCoord4dv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_TexCoord4f(PyObject * self, PyObject * args)
{
    GLfloat arg1, arg2, arg3, arg4, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "ffff", &arg1, &arg2, &arg3, &arg4))
	glTexCoord4f(arg1, arg2, arg3, arg4);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsFloatArray(&op, &vert, &n));
	if (n < 4) {
	    PyErr_SetString(py_gl_Error, "need element with at least 4 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glTexCoord4fv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_TexCoord4i(PyObject * self, PyObject * args)
{
    GLint arg1, arg2, arg3, arg4, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "iiii", &arg1, &arg2, &arg3, &arg4))
	glTexCoord4i(arg1, arg2, arg3, arg4);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsIntArray(&op, &vert, &n));
	if (n < 4) {
	    PyErr_SetString(py_gl_Error, "need element with at least 4 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glTexCoord4iv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_TexCoord4s(PyObject * self, PyObject * args)
{
    GLshort arg1, arg2, arg3, arg4, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "hh", &arg1, &arg2, &arg3, &arg4))
	glTexCoord4s(arg1, arg2, arg3, arg4);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsShortArray(&op, &vert, &n));
	if (n < 4) {
	    PyErr_SetString(py_gl_Error, "need element with at least 4 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glTexCoord4sv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_TexEnvf(PyObject * self, PyObject * args)
{
    GLenum arg1, arg2;
    GLfloat arg3, *rest;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "iif", &arg1, &arg2, &arg3))
	glTexEnvf(arg1, arg2, arg3);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "iiO", &arg1, &arg2, &op));
	TRY(PyArray_AsFloatArray(&op, &rest, &n));
	if ((n == 1) || (n > 3)) {
	    glTexEnvfv(arg1, arg2, rest);
	    PyArray_ClearMemory(op, rest);
	} else {
	    PyErr_SetString(py_gl_Error,
			    "3. argument needs 1 or at least 4 items!");
	    PyArray_ClearMemory(op, rest);
	    return NULL;
	}
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_TexEnvi(PyObject * self, PyObject * args)
{
    GLenum arg1, arg2;
    GLint arg3, *rest;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "iii", &arg1, &arg2, &arg3))
	glTexEnvi(arg1, arg2, arg3);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "iiO", &arg1, &arg2, &op));
	TRY(PyArray_AsIntArray(&op, &rest, &n));
	if ((n == 1) || (n > 3)) {
	    glTexEnviv(arg1, arg2, rest);
	    PyArray_ClearMemory(op, rest);
	} else {
	    PyErr_SetString(py_gl_Error,
			    "3. argument needs 1 or at least 4 items!");
	    PyArray_ClearMemory(op, rest);
	    return NULL;
	}
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_TexGend(PyObject * self, PyObject * args)
{
    GLenum arg1, arg2;
    GLdouble arg3, *rest;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "iid", &arg1, &arg2, &arg3))
	glTexGend(arg1, arg2, arg3);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "iiO", &arg1, &arg2, &op));
	TRY(PyArray_AsDoubleArray(&op, &rest, &n));
	if ((n == 1) || (n > 3)) {
	    glTexGendv(arg1, arg2, rest);
	    PyArray_ClearMemory(op, rest);
	} else {
	    PyErr_SetString(py_gl_Error,
			    "3. argument needs 1 or at least 4 items!");
	    PyArray_ClearMemory(op, rest);
	    return NULL;
	}
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_TexGenf(PyObject * self, PyObject * args)
{
    GLenum arg1, arg2;
    GLfloat arg3, *rest;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "iif", &arg1, &arg2, &arg3))
	glTexGenf(arg1, arg2, arg3);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "iiO", &arg1, &arg2, &op));
	TRY(PyArray_AsFloatArray(&op, &rest, &n));
	if ((n == 1) || (n > 3)) {
	    glTexGenfv(arg1, arg2, rest);
	    PyArray_ClearMemory(op, rest);
	} else {
	    PyErr_SetString(py_gl_Error,
			    "3. argument needs 1 or at least 4 items!");
	    PyArray_ClearMemory(op, rest);
	    return NULL;
	}
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_TexGeni(PyObject * self, PyObject * args)
{
    GLenum arg1, arg2;
    GLint arg3, *rest;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "iii", &arg1, &arg2, &arg3))
	glTexGeni(arg1, arg2, arg3);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "iiO", &arg1, &arg2, &op));
	TRY(PyArray_AsIntArray(&op, &rest, &n));
	if ((n == 1) || (n > 3)) {
	    glTexGeniv(arg1, arg2, rest);
	    PyArray_ClearMemory(op, rest);
	} else {
	    PyErr_SetString(py_gl_Error,
			    "3. argument needs 1 or at least 4 items!");
	    PyArray_ClearMemory(op, rest);
	    return NULL;
	}
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_TexParameterf(PyObject * self, PyObject * args)
{
    GLenum arg1, arg2;
    GLfloat arg3, *rest;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "iif", &arg1, &arg2, &arg3))
	glTexParameterf(arg1, arg2, arg3);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "iiO", &arg1, &arg2, &op));
	TRY(PyArray_AsFloatArray(&op, &rest, &n));
	if ((n == 1) || (n > 3)) {
	    glTexParameterfv(arg1, arg2, rest);
	    PyArray_ClearMemory(op, rest);
	} else {
	    PyErr_SetString(py_gl_Error,
			    "3. argument needs 1 or at least 4 items!");
	    PyArray_ClearMemory(op, rest);
	    return NULL;
	}
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_TexParameteri(PyObject * self, PyObject * args)
{
    GLenum arg1, arg2;
    GLint arg3, *rest;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "iii", &arg1, &arg2, &arg3))
	glTexParameteri(arg1, arg2, arg3);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "iiO", &arg1, &arg2, &op));
	TRY(PyArray_AsIntArray(&op, &rest, &n));
	if ((n == 1) || (n > 3)) {
	    glTexParameteriv(arg1, arg2, rest);
	    PyArray_ClearMemory(op, rest);
	} else {
	    PyErr_SetString(py_gl_Error,
			    "3. argument needs 1 or at least 4 items!");
	    PyArray_ClearMemory(op, rest);
	    return NULL;
	}
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Translated(PyObject * self, PyObject * args)
{
    GLdouble arg1, arg2, arg3;
    TRY(PyArg_ParseTuple(args, "ddd", &arg1, &arg2, &arg3));
    glTranslated(arg1, arg2, arg3);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Translatef(PyObject * self, PyObject * args)
{
    GLfloat arg1, arg2, arg3;
    TRY(PyArg_ParseTuple(args, "fff", &arg1, &arg2, &arg3));
    glTranslatef(arg1, arg2, arg3);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Vertex2d(PyObject * self, PyObject * args)
{
    GLdouble arg1, arg2, *vert;
    PyObject *op;
    int n;
    if (PyArg_ParseTuple(args, "dd", &arg1, &arg2))
	glVertex2d(arg1, arg2);
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsDoubleArray(&op, &vert, &n));
	if (n < 2) {
	    PyErr_SetString(py_gl_Error, "need element with at least 2 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glVertex2dv(vert);
	PyArray_ClearMemory(op, vert);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

/* experimental */

static PyObject *py_gl_ColorVertex2d(PyObject * self, PyObject * args)
{
    GLdouble arg1, arg2, *vert;
    PyObject *op1;
    int n;
#ifdef NUMERIC
    PyObject *op2, *colors = NULL;
    PyArrayObject *cp = NULL;
    int ok1, ok2, nc, m, i;
    GLdouble *x, *y, *c;
#endif
    if (PyArg_ParseTuple(args, "dd", &arg1, &arg2))
	glVertex2d(arg1, arg2);
    else {
	PyErr_Clear();
	if (PyArg_ParseTuple(args, "O", &op1)) {
	    TRY(PyArray_AsDoubleArray(&op1, &vert, &n));
	    if (n < 2) {
		PyErr_SetString(py_gl_Error,
				"need element with at least 2 items");
		PyArray_ClearMemory(op1, vert);
		return NULL;
	    }
	    glVertex2dv(vert);
	    PyArray_ClearMemory(op1, vert);
	}
#ifdef NUMERIC
	else {
	    PyErr_Clear();
	    TRY(PyArg_ParseTuple(args, "OO|O", &op1, &op2, &colors));
	    ok1 = PyArray_AsDoubleArray(&op1, &x, &m);
	    ok2 = PyArray_AsDoubleArray(&op2, &y, &n);
	    if ((m != n) || (!ok1) || (!ok2)) {
		PyErr_SetString(py_gl_Error,
				"coordinate arrays must be of same length or not enough memory");
		PyArray_ClearMemory(op1, x);
		PyArray_ClearMemory(op2, y);
		return NULL;
	    }
	    if (colors) {
		if (!(cp = (PyArrayObject *)
		      PyArray_ContiguousFromObject(colors, PyArray_DOUBLE,
						   1, 2))) {
		    PyArray_ClearMemory(op1, x);
		    PyArray_ClearMemory(op2, y);
		    return NULL;
		}
		c = (GLdouble *) cp->data;
		nc = PyArray_Size((PyObject *) cp);
		if (((nc % 3) != 0) || (n != nc / 3)) {
		    PyErr_SetString(py_gl_Error, "wrong color matrix size");
		    PyArray_ClearMemory(op1, x);
		    PyArray_ClearMemory(op2, y);
		    PyArray_ClearMemory(cp, c);
		    return NULL;
		}
		for (i = 0; i < n; i++) {
		    glColor3dv(c);
		    c += 3;
		    glVertex2d(*x++, *y++);
		}
	    } else {
		for (i = 0; i < n; i++) {
		    glVertex2d(*x++, *y++);
		}
	    }
	    PyArray_ClearMemory(op1, x);
	    PyArray_ClearMemory(op2, y);
	    PyArray_ClearMemory(cp, c);
	}
    }
#else
    }
#endif
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Vertex2f(PyObject * self, PyObject * args)
{
    GLfloat arg1, arg2;
#ifdef NUMERIC
    GLfloat *vert;
    PyObject *op;
    int n;
#endif
    if (PyArg_ParseTuple(args, "ff", &arg1, &arg2))
	glVertex2f(arg1, arg2);
#ifdef NUMERIC
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsFloatArray(&op, &vert, &n));
	if (n < 2) {
	    PyErr_SetString(py_gl_Error, "need element with at least 2 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glVertex2fv(vert);
	PyArray_ClearMemory(op, vert);
    }
#endif
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Vertex2i(PyObject * self, PyObject * args)
{
    GLint arg1, arg2;
#ifdef NUMERIC
    GLint *vert;
    PyObject *op;
    int n;
#endif
    if (PyArg_ParseTuple(args, "ii", &arg1, &arg2))
	glVertex2i(arg1, arg2);
#ifdef NUMERIC
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsIntArray(&op, &vert, &n));
	if (n < 2) {
	    PyErr_SetString(py_gl_Error, "need element with at least 2 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glVertex2iv(vert);
	PyArray_ClearMemory(op, vert);
    }
#endif
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Vertex2s(PyObject * self, PyObject * args)
{
    GLshort arg1, arg2;
#ifdef NUMERIC
    GLshort *vert;
    PyObject *op;
    int n;
#endif
    if (PyArg_ParseTuple(args, "hh", &arg1, &arg2))
	glVertex2s(arg1, arg2);
#ifdef NUMERIC
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsShortArray(&op, &vert, &n));
	if (n < 2) {
	    PyErr_SetString(py_gl_Error, "need element with at least 2 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glVertex2sv(vert);
	PyArray_ClearMemory(op, vert);
    }
#endif
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Vertex3d(PyObject * self, PyObject * args)
{
    GLdouble arg1, arg2, arg3;
#ifdef NUMERIC
    GLdouble *vert;
    PyObject *op;
    int n;
#endif
    if (PyArg_ParseTuple(args, "ddd", &arg1, &arg2, &arg3))
	glVertex3d(arg1, arg2, arg3);
#ifdef NUMERIC
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsDoubleArray(&op, &vert, &n));
	if (n < 3) {
	    PyErr_SetString(py_gl_Error, "need element with at least 3 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glVertex3dv(vert);
	PyArray_ClearMemory(op, vert);
    }
#endif
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Vertex3f(PyObject * self, PyObject * args)
{
    GLfloat arg1, arg2, arg3;
#ifdef NUMERIC
    GLfloat *vert;
    PyObject *op;
    int n;
#endif
    if (PyArg_ParseTuple(args, "fff", &arg1, &arg2, &arg3))
	glVertex3f(arg1, arg2, arg3);
#ifdef NUMERIC
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsFloatArray(&op, &vert, &n));
	if (n < 3) {
	    PyErr_SetString(py_gl_Error, "need element with at least 3 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glVertex3fv(vert);
	PyArray_ClearMemory(op, vert);
    }
#endif
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Vertex3i(PyObject * self, PyObject * args)
{
    GLint arg1, arg2, arg3;
#ifdef NUMERIC
    GLint *vert;
    PyObject *op;
    int n;
#endif
    if (PyArg_ParseTuple(args, "iii", &arg1, &arg2, &arg3))
	glVertex3i(arg1, arg2, arg3);
#ifdef NUMERIC
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsIntArray(&op, &vert, &n));
	if (n < 3) {
	    PyErr_SetString(py_gl_Error, "need element with at least 3 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glVertex3iv(vert);
	PyArray_ClearMemory(op, vert);
    }
#endif
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Vertex3s(PyObject * self, PyObject * args)
{
    GLshort arg1, arg2, arg3;
#ifdef NUMERIC
    GLshort *vert;
    PyObject *op;
    int n;
#endif
    if (PyArg_ParseTuple(args, "hhh", &arg1, &arg2, &arg3))
	glVertex3s(arg1, arg2, arg3);
#ifdef NUMERIC
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsShortArray(&op, &vert, &n));
	if (n < 3) {
	    PyErr_SetString(py_gl_Error, "need element with at least 3 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glVertex3sv(vert);
	PyArray_ClearMemory(op, vert);
    }
#endif
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Vertex4d(PyObject * self, PyObject * args)
{
    GLdouble arg1, arg2, arg3, arg4;
#ifdef NUMERIC
    GLdouble *vert;
    PyObject *op;
    int n;
#endif
    if (PyArg_ParseTuple(args, "dddd", &arg1, &arg2, &arg3, &arg4))
	glVertex4d(arg1, arg2, arg3, arg4);
#ifdef NUMERIC
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsDoubleArray(&op, &vert, &n));
	if (n < 4) {
	    PyErr_SetString(py_gl_Error, "need element with at least 4 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glVertex4dv(vert);
	PyArray_ClearMemory(op, vert);
    }
#endif
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Vertex4f(PyObject * self, PyObject * args)
{
    GLfloat arg1, arg2, arg3, arg4;
#ifdef NUMERIC
    GLfloat *vert;
    PyObject *op;
    int n;
#endif
    if (PyArg_ParseTuple(args, "ffff", &arg1, &arg2, &arg3, &arg4))
	glVertex4f(arg1, arg2, arg3, arg4);
#ifdef NUMERIC
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsFloatArray(&op, &vert, &n));
	if (n < 4) {
	    PyErr_SetString(py_gl_Error, "need element with at least 4 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glVertex4fv(vert);
	PyArray_ClearMemory(op, vert);
    }
#endif
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Vertex4i(PyObject * self, PyObject * args)
{
    GLint arg1, arg2, arg3, arg4;
#ifdef NUMERIC
    GLint *vert;
    PyObject *op;
    int n;
#endif
    if (PyArg_ParseTuple(args, "iiii", &arg1, &arg2, &arg3, &arg4))
	glVertex4i(arg1, arg2, arg3, arg4);
#ifdef NUMERIC
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsIntArray(&op, &vert, &n));
	if (n < 4) {
	    PyErr_SetString(py_gl_Error, "need element with at least 4 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glVertex4iv(vert);
	PyArray_ClearMemory(op, vert);
    }
#endif
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Vertex4s(PyObject * self, PyObject * args)
{
    GLshort arg1, arg2, arg3, arg4;
#ifdef NUMERIC
    GLshort *vert;
    PyObject *op;
    int n;
#endif
    if (PyArg_ParseTuple(args, "hhhh", &arg1, &arg2, &arg3, &arg4))
	glVertex4s(arg1, arg2, arg3, arg4);
#ifdef NUMERIC
    else {
	PyErr_Clear();
	TRY(PyArg_ParseTuple(args, "O", &op));
	TRY(PyArray_AsShortArray(&op, &vert, &n));
	if (n < 4) {
	    PyErr_SetString(py_gl_Error, "need element with at least 4 items");
	    PyArray_ClearMemory(op, vert);
	    return NULL;
	}
	glVertex4sv(vert);
	PyArray_ClearMemory(op, vert);
    }
#endif
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_Viewport(PyObject * self, PyObject * args)
{
    GLint arg1, arg2;
    GLsizei arg3, arg4;
    TRY(PyArg_ParseTuple(args, "iiii", &arg1, &arg2, &arg3, &arg4));
    glViewport(arg1, arg2, arg3, arg4);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_ClipPlane(PyObject * self, PyObject * args)
{
    GLenum arg1;
    PyObject *eop;
    GLdouble *equation;
    int n;
    TRY(PyArg_ParseTuple(args, "iO", &arg1, &eop));
    TRY(PyArray_AsDoubleArray(&eop, &equation, &n));
    if (n < 4) {
	PyErr_SetString(py_gl_Error, "second argument needs at least 4 items");
	Py_DECREF(eop);
	return NULL;
    }
    glClipPlane(arg1, equation);
    Py_DECREF(eop);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_GetClipPlane(PyObject * self, PyObject * args)
{
    GLenum arg1;
#ifdef NUMERIC
    PyArrayObject *values;
    int dims[1];
#else				/* Not NUMERIC */
    PyObject *values;
    int i;
#endif				/* Not NUMERIC */
    GLdouble equation[4];
    TRY(PyArg_ParseTuple(args, "i", &arg1));
    glGetClipPlane(arg1, equation);
#ifdef NUMERIC
    dims[0] = 4;
    values = (PyArrayObject *) PyArray_FromDims(1, dims, PyArray_DOUBLE);
    memcpy(values->data, equation, 4 * sizeof(double));
#else				/* Not NUMERIC */
    values = PyTuple_New(4);
    for (i = 0; i < 4; i++) {
	PyTuple_SET_ITEM(values, i, PyFloat_FromDouble(equation[i]));
    }
#endif				/* Not NUMERIC */
    return (PyObject *) values;
}


static PyObject *py_gl_IsEnabled(PyObject * self, PyObject * args)
{
    GLenum arg1;
    TRY(PyArg_ParseTuple(args, "i", &arg1));
    return PyInt_FromLong(glIsEnabled(arg1));
}

static PyObject *py_gl_IsList(PyObject * self, PyObject * args)
{
    GLuint arg1;
    TRY(PyArg_ParseTuple(args, "i", &arg1));
    return PyInt_FromLong(glIsList(arg1));
}

static PyObject *py_gl_RenderMode(PyObject * self, PyObject * args)
{
    GLenum arg1;
    TRY(PyArg_ParseTuple(args, "i", &arg1));
    return PyInt_FromLong(glRenderMode(arg1));
}

static PyObject *py_gl_GenLists(PyObject * self, PyObject * args)
{
    GLsizei arg1;
    TRY(PyArg_ParseTuple(args, "i", &arg1));
    return PyInt_FromLong(glGenLists(arg1));
}

static PyObject *py_gl_GetDoublev(PyObject * self, PyObject * args)
{
    GLenum pname;
    GLdouble oneparam;
    int nd = 1, dims[3];
#ifdef NUMERIC
    PyArrayObject *params;
#else				/* NUMERIC */
    int i, nitems;
    GLdouble *items;
    PyObject *params;
#endif				/* NUMERIC */
    if (!PyArg_ParseTuple(args, "i", &pname))
	return NULL;
    switch (pname) {
    case GL_ACCUM_CLEAR_VALUE:
#ifdef GL_BLEND_COLOR_EXT
    case GL_BLEND_COLOR_EXT:
#endif
    case GL_COLOR_CLEAR_VALUE:
    case GL_COLOR_WRITEMASK:
    case GL_CURRENT_COLOR:
    case GL_CURRENT_RASTER_COLOR:
    case GL_CURRENT_RASTER_POSITION:
    case GL_CURRENT_RASTER_TEXTURE_COORDS:
    case GL_CURRENT_TEXTURE_COORDS:
    case GL_FOG_COLOR:
    case GL_LIGHT_MODEL_AMBIENT:
    case GL_MAP2_GRID_DOMAIN:
    case GL_SCISSOR_BOX:
    case GL_TEXTURE_ENV_COLOR:
    case GL_VIEWPORT:
	dims[0] = 4;
	break;
    case GL_CURRENT_NORMAL:
	dims[0] = 3;
	break;
    case GL_DEPTH_RANGE:
    case GL_LINE_WIDTH_RANGE:
    case GL_MAP1_GRID_DOMAIN:
    case GL_MAP2_GRID_SEGMENTS:
    case GL_MAX_VIEWPORT_DIMS:
    case GL_POINT_SIZE_RANGE:
    case GL_POLYGON_MODE:
	dims[0] = 2;
	break;
    case GL_MODELVIEW_MATRIX:
    case GL_PROJECTION_MATRIX:
    case GL_TEXTURE_MATRIX:
	nd = 2;
	dims[0] = 4;
	dims[1] = 4;
	break;
    case GL_POLYGON_STIPPLE:
	dims[0] = 32;
	break;
    default:
	glGetDoublev(pname, &oneparam);
	return Py_BuildValue("d", oneparam);
    }
#ifdef NUMERIC
    TRY(params = (PyArrayObject *) PyArray_FromDims(nd, dims, PyArray_DOUBLE));
    glGetDoublev(pname, (double *) params->data);
#else				/* NUMERIC */
    nitems = 1;
    for (i = 0; i < nd; i++) {
	nitems *= dims[i];
    }
    TRY(items = PyMem_NEW(GLdouble, nitems));
    glGetDoublev(pname, items);
    TRY(params = PyTuple_New(nitems));
    for (i = 0; i < nitems; i++) {
	PyTuple_SET_ITEM(params, i, PyFloat_FromDouble(items[i]));
    }
    PyMem_DEL(items); /* items from PyMem_NEW */
#endif				/* NUMERIC */
    return (PyObject *) params;
}

/*OpenGL Extensions */

static PyObject *py_gl_ArrayElementEXT(PyObject * self, PyObject * args)
{
#ifdef GL_EXT_array_element
    GLint arg1;
    TRY(PyArg_ParseTuple(args, "i", &arg1));
    glArrayElementEXT(arg1);
    Py_INCREF(Py_None);
    return Py_None;
#else
    PyErr_SetString(PyExc_ValueError, "ArrayElementEXT not implemented");
    return NULL;
#endif
}

/*These need to be fixed to match the above */
#ifdef NOT_DEFINED

static PyObject *py_gl_PolygonOffsetEXT(PyObject * self, PyObject * args)
{
#ifdef GL_EXT_polygon_offset
    GLfloat arg1, arg2;
    TRY(PyArg_ParseTuple(args, "ff", &arg1, &arg2));
    glPolygonOffsetEXT(arg1, arg2);
    Py_INCREF(Py_None);
    return Py_None;
#else
    PyErr_SetString(PyExc_ValueError, "PolygonOffsetEXT not implemented");
    return NULL;
#endif
}

static PyObject *py_gl_DrawArraysEXT(PyObject * self, PyObject * args)
{
    GLenum arg1;
    GLint arg2;
    GLsizei arg3;
    TRY(PyArg_ParseTuple(args, "iii", &arg1, &arg2, &arg3));
    glDrawArraysEXT(arg1, arg2, arg3);
    Py_INCREF(Py_None);
    return Py_None;
}
static PyObject *py_gl_BlendColorEXT(PyObject * self, PyObject * args)
{
    GLclampf arg1, arg2, arg3, arg4;
    TRY(PyArg_ParseTuple(args, "ffff", &arg1, &arg2, &arg3, &arg4));
    glBlendColorEXT(arg1, arg2, arg3, arg4);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_BlendEquationEXT(PyObject * self, PyObject * args)
{
    GLenum arg1;
    TRY(PyArg_ParseTuple(args, "i", &arg1));
    glBlendEquationEXT(arg1);
    Py_INCREF(Py_None);
    return Py_None;
}

#endif

/*Need to find proper size for glReadPixels array */
static int glformat_size(GLenum format)
{
    switch (format) {
    case GL_COLOR_INDEX:
    case GL_RED:
    case GL_GREEN:
    case GL_BLUE:
    case GL_ALPHA:
    case GL_STENCIL_INDEX:
    case GL_DEPTH_COMPONENT:
    case GL_LUMINANCE:
	return 1;

    case GL_LUMINANCE_ALPHA:
	return 2;

    case GL_RGB:
#ifdef GL_BGR_EXT
    case GL_BGR_EXT:
#endif
	return 3;

    case GL_RGBA:
#ifdef GL_BGRA_EXT
    case GL_BGRA_EXT:
#endif
#ifdef GL_ABGR_EXT
    case GL_ABGR_EXT:
#endif
	return 4;
    default:
	return -1;
    }
}

static int gltype_size(GLenum type)
{
    switch (type) {
    case GL_BYTE:
    case GL_UNSIGNED_BYTE:
	return 8;

    case GL_SHORT:
    case GL_UNSIGNED_SHORT:
	return 16;

    case GL_INT:
    case GL_UNSIGNED_INT:
    case GL_FLOAT:
	return 32;

    case GL_BITMAP:
	return 1;

    default:
	return -1;
    }
}

static PyObject *py_gl_DrawPixels(PyObject * self, PyObject * args)
{
    char *data;
    int format, type, width, height, n;
    int format_size, type_size, total_size;

    TRY(PyArg_ParseTuple(args, "iiiis#", &width, &height, &format, &type,
			 &data, &n));
    /*This is a poor substitute for the RIGHT solution */
    format_size = glformat_size(format);
    ASSERT(format_size != -1, "invalid format");
    type_size = gltype_size(type);
    ASSERT(format_size != -1, "invalid type");
    total_size = type_size * format_size * width * height / 8;
    ASSERT(total_size >= n, "data area too small");

    glDrawPixels(width, height, format, type, data);

    Py_INCREF(Py_None);
    return Py_None;
}

#ifndef MS_WIN32
#ifndef __MWERKS__
#define UNIX
#endif
#endif
#ifdef UNIX
#include <unistd.h>
#include <sys/mman.h>
#endif

#ifdef MS_WIN32
char *my_mmap(char *filename)
#ifdef NUMERIC
;
#else
{
    unsigned long map_size = 0;
    char complaint[256];
    HFILE fh = 0;
    OFSTRUCT file_info;
    HANDLE map_handle;
    char *data;

    fh = OpenFile(filename, &file_info, OF_READWRITE);
    if (fh == HFILE_ERROR) {
	sprintf(complaint, "OpenFile failed: %ld", GetLastError());
	PyErr_SetString(py_gl_Error, complaint);
	return NULL;
    }
    map_size = GetFileSize((HANDLE) fh, NULL);
    map_handle = CreateFileMapping((HANDLE) fh, NULL, PAGE_READWRITE,
				   0, map_size, NULL);
    if (map_handle != NULL) {
	data = (char *) MapViewOfFile(map_handle, FILE_MAP_WRITE, 0, 0, 0);
	if (data != NULL) {
	    return data;
	} else {
	    sprintf(complaint, "MapViewOfFile failed: %ld", GetLastError());
	}
    } else {
	sprintf(complaint, "CreateFileMapping failed: %ld", GetLastError());
    }
    PyErr_SetString(py_gl_Error, complaint);
    return (NULL);
}
#endif
#endif				/* MS_WIN32 */

#if 0
static PyObject *py_gl_DrawMappedPixels(PyObject * self, PyObject * args)
{
    char *data;
    char *filename;
    int format, type, width, height;
    unsigned long n;
    int format_size, type_size, total_size;

    TRY(PyArg_ParseTuple(args, "siiiil", &filename,
			 &width, &height, &format, &type, &n));
    /*This is a poor substitute for the RIGHT solution */
    format_size = glformat_size(format);
    ASSERT(format_size != -1, "invalid format");
    type_size = gltype_size(type);
    ASSERT(format_size != -1, "invalid type");
    total_size = type_size * format_size * width * height / 8;
    ASSERT(total_size >= n, "data area too small");
    data = (char *) my_mmap(filename);
    if (data == NULL) {
	return NULL;
    }
    glDrawPixels(width, height, format, type, data);

    Py_INCREF(Py_None);
    return Py_None;
}

#endif

#ifdef NUMERIC
static PyObject *py_gl_TexImage2D(PyObject * self, PyObject * args)
{
    GLenum arg1;
    GLint arg2;
    GLint arg3;
    GLsizei arg4;
    GLsizei arg5;
    GLint arg6;
    GLenum arg7;
    GLenum arg8;
    PyObject *vop;
    PyArrayObject *ap;
    char *data;
    int format, type;
    int format_size, type_size, total_size;

    TRY(PyArg_ParseTuple
	(args, "iiiiiiiiO", &arg1, &arg2, &arg3, &arg4, &arg5, &arg6,
	 &arg7, &arg8, &vop));

    if (PyString_Check(vop)) {
	data = PyString_AsString(vop);
	format_size = glformat_size(arg7);
	ASSERT(format_size != -1, "invalid format");
	type_size = gltype_size(arg8);
	ASSERT(format_size != -1, "invalid type");
	total_size = type_size * format_size * arg4 * arg5 / 8;
	ASSERT(total_size >= PyString_Size(vop), "data area too small");

	glTexImage2D(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, data);
    } else {
	if (PyArray_Check(vop)) {
	    ap = (PyArrayObject *) vop;
	    ASSERT(ap->nd == 2 || ap->nd == 3, "array must be either 2 or 3d");
	} else {
	    TRY(ap =
		(PyArrayObject *) PyArray_ContiguousFromObject(vop,
							       PyArray_UBYTE,
							       2, 3));
	}

	type = typecode2gltype[ap->descr->type_num];
	ASSERT(type != -1, "can't convert this type of array to an image");

	if (ap->nd == 2) {
	    format = GL_LUMINANCE;
	} else {
	    ASSERT(ap->dimensions[2] == 3
		   || ap->dimensions[2] == 4, "3d array must be RGB or RGBA");

	    if (ap->dimensions[2] == 3) {
		format = GL_RGB;
	    } else {
		format = GL_RGBA;
	    }
	}
	glTexImage2D(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, ap->data);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_TexImage1D(PyObject * self, PyObject * args)
{
    GLenum arg1;
    GLint arg2;
    GLint arg3;
    GLsizei arg4;
    GLint arg6;
    GLenum arg7;
    GLenum arg8;
    PyObject *vop;
    PyArrayObject *ap;
    char *data;
    int format, type;
    int format_size, type_size, total_size;

    TRY(PyArg_ParseTuple
	(args, "iiiiiiiO", &arg1, &arg2, &arg3, &arg4, &arg6, &arg7,
	 &arg8, &vop));

    if (PyString_Check(vop)) {
	data = PyString_AsString(vop);
	format_size = glformat_size(arg7);
	ASSERT(format_size != -1, "invalid format");
	type_size = gltype_size(arg8);
	ASSERT(format_size != -1, "invalid type");
	total_size = type_size * format_size * arg4 / 8;
	ASSERT(total_size >= PyString_Size(vop), "data area too small");

	glTexImage1D(arg1, arg2, arg3, arg4, arg6, arg7, arg8, data);
    } else {
	if (PyArray_Check(vop)) {
	    ap = (PyArrayObject *) vop;
	    ASSERT(ap->nd == 2, "array must be either 2d");
	} else {
	    TRY(ap =
		(PyArrayObject *) PyArray_ContiguousFromObject(vop,
							       PyArray_UBYTE,
							       2, 3));
	}

	type = typecode2gltype[ap->descr->type_num];
	ASSERT(type != -1, "can't convert this type of array to an image");

	ASSERT(ap->dimensions[1] == 3
	       || ap->dimensions[1] == 4, "3d array must be RGB or RGBA");

	if (ap->dimensions[1] == 3) {
	    format = GL_RGB;
	} else {
	    format = GL_RGBA;
	}
	glTexImage1D(arg1, arg2, arg3, arg4, arg6, arg7, arg8, ap->data);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

#endif

static PyObject *py_gl_ReadPixels(PyObject * self, PyObject * args)
{
    PyObject *sp;
    char *data;
    int x, y;
    int format, type, width, height, n = -1;
    int format_size, type_size, total_size;

    TRY(PyArg_ParseTuple
	(args, "iiiiii|i", &x, &y, &width, &height, &format, &type, &n));
    /*This is a poor substitute for the RIGHT solution */
    format_size = glformat_size(format);
    ASSERT(format_size != -1, "invalid format");
    type_size = gltype_size(type);
    ASSERT(format_size != -1, "invalid type");
    total_size = type_size * format_size * width * height / 8;
    if (n == -1)
	n = total_size;
    ASSERT(total_size >= n, "data area too small");

    data = (char *) malloc(n);
    glReadPixels(x, y, width, height, format, type, data);
    sp = PyString_FromStringAndSize(data, n);
    free(data);
    return sp;
}

#ifdef GL_VERSION_1_1
#define OPENGL1_1
#endif

/* 
   ############################### OpenGL 1.1 ########################
 */
#ifdef OPENGL1_1
#ifdef NUMERIC

static int type2gl[] =
    { -1, GL_UNSIGNED_BYTE, GL_BYTE, GL_SHORT, GL_INT, -1, GL_FLOAT,
    GL_DOUBLE, -1, -1, -1, -1
};

static PyObject *py_gl_VertexPointer(PyObject * self, PyObject * args)
{
    GLint size;
    GLenum type;
    GLsizei stride;
    PyObject *op;
    PyArrayObject *ap;
    TRY(PyArg_ParseTuple(args, "iiO", &size, &stride, &op));
    if (PyArray_Check(op))
	ap = (PyArrayObject *) op;
    else {
	TRY(ap =
	    (PyArrayObject *) PyArray_ContiguousFromObject(op,
							   PyArray_DOUBLE,
							   1, 0));
    }
    type = type2gl[ap->descr->type_num];
    ASSERT(type != -1, "Can't convert this type of array!");
    glVertexPointer(size, type, stride, ap->data);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_ColorPointer(PyObject * self, PyObject * args)
{
    GLint size;
    GLenum type;
    GLsizei stride;
    PyObject *op;
    PyArrayObject *ap;
    TRY(PyArg_ParseTuple(args, "iiO", &size, &stride, &op));
    if (PyArray_Check(op))
	ap = (PyArrayObject *) op;
    else {
	TRY(ap =
	    (PyArrayObject *) PyArray_ContiguousFromObject(op,
							   PyArray_DOUBLE,
							   1, 0));
    }
    type = type2gl[ap->descr->type_num];
    ASSERT(type != -1, "Can't convert this type of array!");
    glColorPointer(size, type, stride, ap->data);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_NormalPointer(PyObject * self, PyObject * args)
{
    GLint size;
    GLenum type;
    GLsizei stride;
    PyObject *op;
    PyArrayObject *ap;
    TRY(PyArg_ParseTuple(args, "iiO", &size, &stride, &op));
    if (PyArray_Check(op))
	ap = (PyArrayObject *) op;
    else {
	TRY(ap =
	    (PyArrayObject *) PyArray_ContiguousFromObject(op,
							   PyArray_DOUBLE,
							   1, 0));
    }
    type = type2gl[ap->descr->type_num];
    ASSERT(type != -1, "Can't convert this type of array!");
    glNormalPointer(type, stride, ap->data);
    Py_INCREF(Py_None);
    return Py_None;
}
static PyObject *py_gl_TexCoordPointer(PyObject * self, PyObject * args)
{
    GLint size;
    GLenum type;
    GLsizei stride;
    PyObject *op;
    PyArrayObject *ap;
    TRY(PyArg_ParseTuple(args, "iiO", &size, &stride, &op));
    if (PyArray_Check(op))
	ap = (PyArrayObject *) op;
    else {
	TRY(ap =
	    (PyArrayObject *) PyArray_ContiguousFromObject(op,
							   PyArray_DOUBLE,
							   1, 0));
    }
    type = type2gl[ap->descr->type_num];
    ASSERT(type != -1, "Can't convert this type of array!");
    glTexCoordPointer(size, type, stride, ap->data);
    Py_INCREF(Py_None);
    return Py_None;
}
static PyObject *py_gl_DrawArrays(PyObject * self, PyObject * args)
{
    GLenum mode;
    GLint first;
    GLsizei count;
    TRY(PyArg_ParseTuple(args, "iii", &mode, &first, &count));
    glDrawArrays(mode, first, count);
    Py_INCREF(Py_None);
    return Py_None;
}

#endif				/* NUMERIC */
#endif				/* OPENGL1_1 */

static PyObject *py_gl_GenTextures(PyObject * self, PyObject * args)
{
    GLsizei arg1;
    GLuint texindex;
    TRY(PyArg_ParseTuple(args, "i", &arg1));
    glGenTextures(arg1, &texindex);
    return PyInt_FromLong(texindex);
}

static PyObject *py_gl_BindTexture(PyObject * self, PyObject * args)
{
    GLenum arg1;
    GLint arg2;
    TRY(PyArg_ParseTuple(args, "ii", &arg1, &arg2));
    glBindTexture(arg1, arg2);
    Py_INCREF(Py_None);
    return Py_None;
}

/* Michel Sanner Jan 14 1998 */

/*
   Check for GL_ERRORS and returns error code or None
 */
static PyObject *py_gl_GetError(PyObject * self, PyObject * args)
{
    GLenum errCode = glGetError();

    if (errCode == GL_NO_ERROR) {
	Py_INCREF(Py_None);
	return Py_None;
    } else {
	return Py_BuildValue("i", errCode);
    }
}

#ifdef NUMERIC
/*
   Allocate a numeric array of a user specified size to hold hardware picking
   event hits. This "Numeric" array is return to the interpreter.
   OpenGL 1.0, requires Numeric extension.

   available from the interpreter as:
   SelectBuffer( int size )
 */
static PyObject *py_gl_SelectBuffer(PyObject * self, PyObject * args)
{
    int size;
    PyArrayObject *buf;

    TRY(PyArg_ParseTuple(args, "i", &size));
    TRY(buf = (PyArrayObject *) PyArray_FromDims(1, &size, PyArray_INT));

    glSelectBuffer(buf->dimensions[0], (GLuint *) buf->data);

    return (PyObject *) buf;
}

#endif


/* END Michel Sanner */

/* calling glBitmap with only 6 arguments means calling it with NULL as the seventh -- see Ref Manual, 1.1, page 313 */

static PyObject *py_gl_Bitmap(PyObject * self, PyObject * args)
{
    int num_args;
    int width, height;
    float xorig, yorig, xmove, ymove;
    unsigned char *bitmap;
    int length;

    num_args = PyArg_ParseTuple(args, "iiffff|s#",
				&width, &height,
				&xorig, &yorig,
				&xmove, &ymove, &bitmap, &length);
    switch (num_args) {
    case 0:
	return NULL;
    case 6:
	bitmap = NULL;		/* no break on purpose */
    default:
	glBitmap(width, height, xorig, yorig, xmove, ymove, bitmap);
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_gl_LoadMatrix(PyObject * self, PyObject * args)
{
    PyObject *op;
    GLdouble *vert;
    int n;
    TRY(PyArg_ParseTuple(args, "O", &op));
    TRY(PyArray_AsDoubleArray(&op, &vert, &n));
    if (n != 16) {
	PyErr_SetString(py_gl_Error, "argument should have 16 values");
	PyArray_ClearMemory(op, vert);
	return NULL;
    }
    glLoadMatrixd(vert);
    Py_INCREF(Py_None);
    return Py_None;
}

/* 
   ########################################################################
 */

static PyMethodDef py_gl_methods[] = {
    {"glAccum", py_gl_Accum, 1},
    {"glAlphaFunc", py_gl_AlphaFunc, 1},
    {"glBegin", py_gl_Begin, 1},
    {"glBlendFunc", py_gl_BlendFunc, 1},
    {"glCallList", py_gl_CallList, 1},
    {"glClear", py_gl_Clear, 1},
    {"glClearAccum", py_gl_ClearAccum, 1},
    {"glClearColor", py_gl_ClearColor, 1},
    {"glClearDepth", py_gl_ClearDepth, 1},
    {"glClearIndex", py_gl_ClearIndex, 1},
    {"glClearStencil", py_gl_ClearStencil, 1},
    {"glColor3", py_gl_Color3d, 1},
    {"glColor3b", py_gl_Color3b, 1},
    {"glColor3bv", py_gl_Color3b, 1},
    {"glColor3d", py_gl_Color3d, 1},
    {"glColor3dv", py_gl_Color3d, 1},
    {"glColor3f", py_gl_Color3f, 1},
    {"glColor3fv", py_gl_Color3f, 1},
    {"glColor3i", py_gl_Color3i, 1},
    {"glColor3iv", py_gl_Color3i, 1},
    {"glColor3s", py_gl_Color3s, 1},
    {"glColor3sv", py_gl_Color3s, 1},
    {"glColor3ub", py_gl_Color3ub, 1},
    {"glColor3ubv", py_gl_Color3ub, 1},
    {"glColor3ui", py_gl_Color3ui, 1},
    {"glColor3uiv", py_gl_Color3ui, 1},
    {"glColor3us", py_gl_Color3us, 1},
    {"glColor3usv", py_gl_Color3us, 1},
    {"glColor4", py_gl_Color4d, 1},
    {"glColor4b", py_gl_Color4b, 1},
    {"glColor4bv", py_gl_Color4b, 1},
    {"glColor4d", py_gl_Color4d, 1},
    {"glColor4dv", py_gl_Color4d, 1},
    {"glColor4f", py_gl_Color4f, 1},
    {"glColor4fv", py_gl_Color4f, 1},
    {"glColor4i", py_gl_Color4i, 1},
    {"glColor4iv", py_gl_Color4i, 1},
    {"glColor4s", py_gl_Color4s, 1},
    {"glColor4sv", py_gl_Color4s, 1},
    {"glColor4ub", py_gl_Color4ub, 1},
    {"glColor4ubv", py_gl_Color4ub, 1},
    {"glColor4ui", py_gl_Color4ui, 1},
    {"glColor4uiv", py_gl_Color4ui, 1},
    {"glColor4us", py_gl_Color4us, 1},
    {"glColor4usv", py_gl_Color4us, 1},
    {"glColorMask", py_gl_ColorMask, 1},
    {"glColorMaterial", py_gl_ColorMaterial, 1},
    {"glCopyPixels", py_gl_CopyPixels, 1},
    {"glCullFace", py_gl_CullFace, 1},
    {"glDeleteLists", py_gl_DeleteLists, 1},
    {"glDepthFunc", py_gl_DepthFunc, 1},
    {"glDepthMask", py_gl_DepthMask, 1},
    {"glDepthRange", py_gl_DepthRange, 1},
    {"glDisable", py_gl_Disable, 1},
    {"glDrawBuffer", py_gl_DrawBuffer, 1},
    {"glEdgeFlag", py_gl_EdgeFlag, 1},
    {"glEnable", py_gl_Enable, 1},
    {"glEnd", py_gl_End, 1},
    {"glEndList", py_gl_EndList, 1},
    {"glEvalCoord1", py_gl_EvalCoord1d, 1},
    {"glEvalCoord1d", py_gl_EvalCoord1d, 1},
    {"glEvalCoord1dv", py_gl_EvalCoord1d, 1},
    {"glEvalCoord1f", py_gl_EvalCoord1f, 1},
    {"glEvalCoord1fv", py_gl_EvalCoord1f, 1},
    {"glEvalCoord2", py_gl_EvalCoord2d, 1},
    {"glEvalCoord2d", py_gl_EvalCoord2d, 1},
    {"glEvalCoord2dv", py_gl_EvalCoord2d, 1},
    {"glEvalCoord2f", py_gl_EvalCoord2f, 1},
    {"glEvalCoord2fv", py_gl_EvalCoord2f, 1},
    {"glEvalMesh1", py_gl_EvalMesh1, 1},
    {"glEvalMesh2", py_gl_EvalMesh2, 1},
    {"glEvalPoint1", py_gl_EvalPoint1, 1},
    {"glEvalPoint2", py_gl_EvalPoint2, 1},
    {"glFinish", py_gl_Finish, 1},
    {"glFlush", py_gl_Flush, 1},
    {"glFog", py_gl_Fogf, 1},
    {"glFogf", py_gl_Fogf, 1},
    {"glFogfv", py_gl_Fogf, 1},
    {"glFogi", py_gl_Fogi, 1},
    {"glFogiv", py_gl_Fogi, 1},
    {"glFrontFace", py_gl_FrontFace, 1},
    {"glFrustum", py_gl_Frustum, 1},
    {"glGetString", py_gl_GetString, 1},
    {"glHint", py_gl_Hint, 1},
    {"glIndexMask", py_gl_IndexMask, 1},
    {"glIndex", py_gl_Indexd, 1},
    {"glIndexd", py_gl_Indexd, 1},
    {"glIndexdv", py_gl_Indexd, 1},
    {"glIndexf", py_gl_Indexf, 1},
    {"glIndexfv", py_gl_Indexf, 1},
    {"glIndexi", py_gl_Indexi, 1},
    {"glIndexiv", py_gl_Indexi, 1},
    {"glIndexs", py_gl_Indexs, 1},
    {"glIndexsv", py_gl_Indexs, 1},
    {"glInitNames", py_gl_InitNames, 1},
    {"glLightModel", py_gl_LightModelf, 1},
    {"glLightModelf", py_gl_LightModelf, 1},
    {"glLightModelfv", py_gl_LightModelf, 1},
    {"glLightModeli", py_gl_LightModeli, 1},
    {"glLightModeliv", py_gl_LightModeli, 1},
    {"glLight", py_gl_Lightf, 1},
    {"glLightf", py_gl_Lightf, 1},
    {"glLightfv", py_gl_Lightf, 1},
    {"glLighti", py_gl_Lighti, 1},
    {"glLightiv", py_gl_Lighti, 1},
    {"glLineStipple", py_gl_LineStipple, 1},
    {"glLineWidth", py_gl_LineWidth, 1},
    {"glListBase", py_gl_ListBase, 1},
    {"glLoadIdentity", py_gl_LoadIdentity, 1},
    {"glLoadName", py_gl_LoadName, 1},
    {"glLogicOp", py_gl_LogicOp, 1},
    {"glMap1d", py_gl_Map1d, 1},
    {"glMap1f", py_gl_Map1f, 1},
    {"glMap2d", py_gl_Map2d, 1},
    {"glMap2f", py_gl_Map2f, 1},
    {"glMapGrid1d", py_gl_MapGrid1d, 1},
    {"glMapGrid1f", py_gl_MapGrid1f, 1},
    {"glMapGrid2d", py_gl_MapGrid2d, 1},
    {"glMapGrid2f", py_gl_MapGrid2f, 1},
    {"glMaterial", py_gl_Materialf, 1},
    {"glMaterialf", py_gl_Materialf, 1},
    {"glMaterialfv", py_gl_Materialf, 1},
    {"glMateriali", py_gl_Materiali, 1},
    {"glMaterialiv", py_gl_Materiali, 1},
    {"glMatrixMode", py_gl_MatrixMode, 1},
    {"glMultMatrix", py_gl_MultMatrixd, 1},
    {"glMultMatrixd", py_gl_MultMatrixd, 1},
    {"glMultMatrixf", py_gl_MultMatrixf, 1},
    {"glNewList", py_gl_NewList, 1},
    {"glNormal3", py_gl_Normal3d, 1},
    {"glNormal3b", py_gl_Normal3s, 1},
    {"glNormal3bv", py_gl_Normal3s, 1},
    {"glNormal3d", py_gl_Normal3d, 1},
    {"glNormal3dv", py_gl_Normal3d, 1},
    {"glNormal3f", py_gl_Normal3f, 1},
    {"glNormal3fv", py_gl_Normal3f, 1},
    {"glNormal3i", py_gl_Normal3i, 1},
    {"glNormal3iv", py_gl_Normal3i, 1},
    {"glNormal3s", py_gl_Normal3s, 1},
    {"glNormal3sv", py_gl_Normal3s, 1},
    {"glOrtho", py_gl_Ortho, 1},
    {"glPassThrough", py_gl_PassThrough, 1},
    {"glPixelStoref", py_gl_PixelStoref, 1},
    {"glPixelStorei", py_gl_PixelStorei, 1},
    {"glPixelTransferf", py_gl_PixelTransferf, 1},
    {"glPixelTransferi", py_gl_PixelTransferi, 1},
    {"glPixelZoom", py_gl_PixelZoom, 1},
    {"glPointSize", py_gl_PointSize, 1},
    {"glPolygonMode", py_gl_PolygonMode, 1},
#ifdef NUMERIC
    {"glPolygonStipple", py_gl_PolygonStipple, 1},
#endif
    {"glPopAttrib", py_gl_PopAttrib, 1},
    {"glPopMatrix", py_gl_PopMatrix, 1},
    {"glPopName", py_gl_PopName, 1},
    {"glPushAttrib", py_gl_PushAttrib, 1},
    {"glPushMatrix", py_gl_PushMatrix, 1},
    {"glPushName", py_gl_PushName, 1},
    {"glRasterPos2", py_gl_RasterPos2d, 1},
    {"glRasterPos2d", py_gl_RasterPos2d, 1},
    {"glRasterPos2dv", py_gl_RasterPos2d, 1},
    {"glRasterPos2f", py_gl_RasterPos2f, 1},
    {"glRasterPos2fv", py_gl_RasterPos2f, 1},
    {"glRasterPos2i", py_gl_RasterPos2i, 1},
    {"glRasterPos2iv", py_gl_RasterPos2i, 1},
    {"glRasterPos2s", py_gl_RasterPos2s, 1},
    {"glRasterPos2sv", py_gl_RasterPos2s, 1},
    {"glRasterPos3", py_gl_RasterPos3d, 1},
    {"glRasterPos3d", py_gl_RasterPos3d, 1},
    {"glRasterPos3dv", py_gl_RasterPos3d, 1},
    {"glRasterPos3f", py_gl_RasterPos3f, 1},
    {"glRasterPos3fv", py_gl_RasterPos3f, 1},
    {"glRasterPos3i", py_gl_RasterPos3i, 1},
    {"glRasterPos3iv", py_gl_RasterPos3i, 1},
    {"glRasterPos3s", py_gl_RasterPos3s, 1},
    {"glRasterPos3sv", py_gl_RasterPos3s, 1},
    {"glRasterPos4", py_gl_RasterPos4d, 1},
    {"glRasterPos4d", py_gl_RasterPos4d, 1},
    {"glRasterPos4dv", py_gl_RasterPos4d, 1},
    {"glRasterPos4f", py_gl_RasterPos4f, 1},
    {"glRasterPos4fv", py_gl_RasterPos4f, 1},
    {"glRasterPos4i", py_gl_RasterPos4i, 1},
    {"glRasterPos4iv", py_gl_RasterPos4i, 1},
    {"glRasterPos4s", py_gl_RasterPos4s, 1},
    {"glRasterPos4sv", py_gl_RasterPos4s, 1},
    {"glReadBuffer", py_gl_ReadBuffer, 1},
    {"glRect", py_gl_Rectd, 1},
    {"glRectd", py_gl_Rectd, 1},
    {"glRectdv", py_gl_Rectd, 1},
    {"glRectf", py_gl_Rectf, 1},
    {"glRectfv", py_gl_Rectf, 1},
    {"glRecti", py_gl_Recti, 1},
    {"glRectiv", py_gl_Recti, 1},
    {"glRects", py_gl_Rects, 1},
    {"glRectsv", py_gl_Rects, 1},
    {"glRotate", py_gl_Rotated, 1},
    {"glRotated", py_gl_Rotated, 1},
    {"glRotated", py_gl_Rotated, 1},
    {"glRotatef", py_gl_Rotatef, 1},
    {"glScale", py_gl_Scaled, 1},
    {"glScaled", py_gl_Scaled, 1},
    {"glScalef", py_gl_Scalef, 1},
    {"glScissor", py_gl_Scissor, 1},
    {"glShadeModel", py_gl_ShadeModel, 1},
    {"glStencilFunc", py_gl_StencilFunc, 1},
    {"glStencilMask", py_gl_StencilMask, 1},
    {"glStencilOp", py_gl_StencilOp, 1},
    {"glTexCoord1", py_gl_TexCoord1d, 1},
    {"glTexCoord1d", py_gl_TexCoord1d, 1},
    {"glTexCoord1dv", py_gl_TexCoord1d, 1},
    {"glTexCoord1f", py_gl_TexCoord1f, 1},
    {"glTexCoord1fv", py_gl_TexCoord1f, 1},
    {"glTexCoord1i", py_gl_TexCoord1i, 1},
    {"glTexCoord1iv", py_gl_TexCoord1i, 1},
    {"glTexCoord1s", py_gl_TexCoord1s, 1},
    {"glTexCoord1sv", py_gl_TexCoord1s, 1},
    {"glTexCoord2", py_gl_TexCoord2d, 1},
    {"glTexCoord2d", py_gl_TexCoord2d, 1},
    {"glTexCoord2dv", py_gl_TexCoord2d, 1},
    {"glTexCoord2f", py_gl_TexCoord2f, 1},
    {"glTexCoord2fv", py_gl_TexCoord2f, 1},
    {"glTexCoord2i", py_gl_TexCoord2i, 1},
    {"glTexCoord2iv", py_gl_TexCoord2i, 1},
    {"glTexCoord2s", py_gl_TexCoord2s, 1},
    {"glTexCoord2sv", py_gl_TexCoord2s, 1},
    {"glTexCoord3", py_gl_TexCoord3d, 1},
    {"glTexCoord3d", py_gl_TexCoord3d, 1},
    {"glTexCoord3dv", py_gl_TexCoord3d, 1},
    {"glTexCoord3f", py_gl_TexCoord3f, 1},
    {"glTexCoord3fv", py_gl_TexCoord3f, 1},
    {"glTexCoord3i", py_gl_TexCoord3i, 1},
    {"glTexCoord3iv", py_gl_TexCoord3i, 1},
    {"glTexCoord3s", py_gl_TexCoord3s, 1},
    {"glTexCoord3sv", py_gl_TexCoord3s, 1},
    {"glTexCoord4", py_gl_TexCoord4d, 1},
    {"glTexCoord4d", py_gl_TexCoord4d, 1},
    {"glTexCoord4dv", py_gl_TexCoord4d, 1},
    {"glTexCoord4f", py_gl_TexCoord4f, 1},
    {"glTexCoord4fv", py_gl_TexCoord4f, 1},
    {"glTexCoord4i", py_gl_TexCoord4i, 1},
    {"glTexCoord4iv", py_gl_TexCoord4i, 1},
    {"glTexCoord4s", py_gl_TexCoord4s, 1},
    {"glTexCoord4sv", py_gl_TexCoord4s, 1},
    {"glTexEnv", py_gl_TexEnvf, 1},
    {"glTexEnvf", py_gl_TexEnvf, 1},
    {"glTexEnvfv", py_gl_TexEnvf, 1},
    {"glTexEnvi", py_gl_TexEnvi, 1},
    {"glTexEnviv", py_gl_TexEnvi, 1},
    {"glTexGen", py_gl_TexGend, 1},
    {"glTexGend", py_gl_TexGend, 1},
    {"glTexGendv", py_gl_TexGend, 1},
    {"glTexGenf", py_gl_TexGenf, 1},
    {"glTexGenfv", py_gl_TexGenf, 1},
    {"glTexGeni", py_gl_TexGeni, 1},
    {"glTexGeniv", py_gl_TexGeni, 1},
    {"glTexParameter", py_gl_TexParameterf, 1},
    {"glTexParameterf", py_gl_TexParameterf, 1},
    {"glTexParameterfv", py_gl_TexParameterf, 1},
    {"glTexParameteri", py_gl_TexParameteri, 1},
    {"glTexParameteriv", py_gl_TexParameteri, 1},
    {"glTranslate", py_gl_Translated, 1},
    {"glTranslated", py_gl_Translated, 1},
    {"glTranslatef", py_gl_Translatef, 1},
    {"glVertex2", py_gl_Vertex2d, 1},
    {"glVertex2d", py_gl_Vertex2d, 1},
    {"glVertex2dv", py_gl_Vertex2d, 1},
    {"glColorVertex2", py_gl_ColorVertex2d, 1},
    {"glColorVertex2d", py_gl_ColorVertex2d, 1},
    {"glVertex2f", py_gl_Vertex2f, 1},
    {"glVertex2fv", py_gl_Vertex2f, 1},
    {"glVertex2i", py_gl_Vertex2i, 1},
    {"glVertex2iv", py_gl_Vertex2i, 1},
    {"glVertex2s", py_gl_Vertex2s, 1},
    {"glVertex2sv", py_gl_Vertex2s, 1},
    {"glVertex3", py_gl_Vertex3d, 1},
    {"glVertex3d", py_gl_Vertex3d, 1},
    {"glVertex3dv", py_gl_Vertex3d, 1},
    {"glVertex3f", py_gl_Vertex3f, 1},
    {"glVertex3fv", py_gl_Vertex3f, 1},
    {"glVertex3i", py_gl_Vertex3i, 1},
    {"glVertex3iv", py_gl_Vertex3i, 1},
    {"glVertex3s", py_gl_Vertex3s, 1},
    {"glVertex3sv", py_gl_Vertex3s, 1},
    {"glVertex4", py_gl_Vertex4d, 1},
    {"glVertex4d", py_gl_Vertex4d, 1},
    {"glVertex4dv", py_gl_Vertex4d, 1},
    {"glVertex4f", py_gl_Vertex4f, 1},
    {"glVertex4fv", py_gl_Vertex4f, 1},
    {"glVertex4i", py_gl_Vertex4i, 1},
    {"glVertex4iv", py_gl_Vertex4i, 1},
    {"glVertex4s", py_gl_Vertex4s, 1},
    {"glVertex4sv", py_gl_Vertex4s, 1},
    {"glViewport", py_gl_Viewport, 1},

    {"glClipPlane", py_gl_ClipPlane, 1},
    {"glGetClipPlane", py_gl_GetClipPlane, 1},
    {"glIsEnabled", py_gl_IsEnabled, 1},
    {"glIsList", py_gl_IsList, 1},
    {"glRenderMode", py_gl_RenderMode, 1},
    {"glGenLists", py_gl_GenLists, 1},
    {"glGetDouble", py_gl_GetDoublev, 1},
    {"glGetDoublev", py_gl_GetDoublev, 1},
    {"glDrawPixels", py_gl_DrawPixels, 1},
#if 0
    {"glDrawMappedPixels", py_gl_DrawMappedPixels, 1},
#endif
#ifdef NUMERIC
    {"glTexImage2D", py_gl_TexImage2D, 1},
    {"glTexImage1D", py_gl_TexImage1D, 1},
    {"glSelectBuffer", py_gl_SelectBuffer, 1},
#endif				/* NUMERIC */
    {"glGetError", py_gl_GetError, 1},
    {"glReadPixels", py_gl_ReadPixels, 1},

/*Various OpenGL extensions */
    {"glArrayElementEXT", py_gl_ArrayElementEXT, 1},
#ifdef NOT_DEFINED
    {"glPolygonOffsetEXT", py_gl_PolygonOffsetEXT, 1},
    {"glDrawArraysEXT", py_gl_DrawArraysEXT, 1},
    {"glBlendColorEXT", py_gl_BlendColorEXT, 1},
    {"glBlendEquationEXT", py_gl_BlendEquationEXT, 1},
#endif

#ifdef OPENGL1_1
    {"glGenTextures", py_gl_GenTextures, 1},
    {"glBindTexture", py_gl_BindTexture, 1},
#ifdef NUMERIC
    {"glVertexPointer", py_gl_VertexPointer, 1},
    {"glColorPointer", py_gl_ColorPointer, 1},
    {"glNormalPointer", py_gl_NormalPointer, 1},
    {"glTexCoordPointer", py_gl_TexCoordPointer, 1},
    {"glDrawArrays", py_gl_DrawArrays, 1},
#endif				/* NUMERIC */
#endif				/* OPENGL1_1 */
    {"glBitmap", py_gl_Bitmap, 1},
    {"glLoadMatrix", py_gl_LoadMatrix, 1},
    {NULL, NULL}
};


#ifdef NUMERIC
DL_EXPORT(void) init_opengl_num(void);
DL_EXPORT(void) init_opengl_num(void)
#else
DL_EXPORT(void) init_opengl(void);
DL_EXPORT(void) init_opengl(void)
#endif
{
    PyObject *m, *d;
    PyObject *nl, *gl;

#ifdef NUMERIC
    m = Py_InitModule("_opengl_num", py_gl_methods);

#ifdef import_array		/* Konrad Hinsen's version */
    import_array();
#endif
#else
    m = Py_InitModule("_opengl", py_gl_methods);
#endif
    d = PyModule_GetDict(m);
#ifdef NUMERIC
    py_gl_Error = Py_BuildValue("s", "_opengl_num.error");
#else
    py_gl_Error = Py_BuildValue("s", "_opengl.error");
#endif
    PyDict_SetItemString(d, "error", py_gl_Error);
#ifdef NUMERIC
    nl = PyInt_FromLong(1L);
#else
    nl = PyInt_FromLong(0L);
#endif
    PyDict_SetItemString(d, "_numeric", nl);
    Py_DECREF(nl);
#ifndef _PYMOL_NO_GLUT
    gl = PyInt_FromLong(1L);
#else
    gl = PyInt_FromLong(0L);
#endif
    PyDict_SetItemString(d, "_glut", gl);
    Py_DECREF(gl);
    if (PyErr_Occurred())
#ifdef NUMERIC
	Py_FatalError("can't initialize module _opengl_num");
#else
	Py_FatalError("can't initialize module _opengl");
#endif
}

/* for distutils compatibility on WIN32 */
#ifndef NUMERIC
void init_openglmodule(void);
void init_openglmodule(void) {init_opengl();}
#endif
