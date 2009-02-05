#ifndef OPENGL_UTIL_H
#define OPENGL_UTIL_H

static PyObject *gl_Error;
#define TRY(E) if(! (E)) return NULL
#define ASSERT(E,M) if(!(E)) return ErrorReturn(M)

#ifdef NUMERIC
static int typecode2gltype[] =
{-1, GL_UNSIGNED_BYTE, GL_BYTE, GL_SHORT, GL_INT, -1, GL_FLOAT, -1, -1, -1, -1, -1};

#if   defined(HAVE_ARRAYOBJECT_H)
#include "arrayobject.h"
#elif defined(HAVE_EXTENSIONS_ARRAYOBJECT_H)
#include "Extensions/arrayobject.h"
#elif defined(HAVE_NUMERIC_ARRAYOBJECT_H)
#include "Numeric/arrayobject.h"
#elif defined(HAVE_NUMERICAL_ARRAYOBJECT_H)
#include "numerical/arrayobject.h"
#else
/*#error "Don't know where to find file 'arrayobject.h'" */
#include "arrayobject.h"
#endif

#define PyArray_ClearMemory(op, pitems) Py_DECREF((op))

static int PyArray_AsDoubleArray(PyObject ** op, GLdouble ** pitems, int *pn)
{
    PyArrayObject *mp;
    mp = (PyArrayObject *) PyArray_ContiguousFromObject(*op, PyArray_DOUBLE, 0, 1);
    if (!mp)
	return 0;
    *pitems = (GLdouble *) (mp->data);
    *pn = mp->dimensions[0];
    *op = (PyObject *) mp;
    return 1;
}

/* not used anywhere
static int PyArray_AsFloatArray(PyObject ** op, GLfloat ** pitems, int *pn)
{
    PyArrayObject *mp;
    mp = (PyArrayObject *) PyArray_ContiguousFromObject(*op, PyArray_FLOAT, 0, 1);
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
    mp = (PyArrayObject *) PyArray_ContiguousFromObject(*op, PyArray_SHORT, 0, 1);
    if (!mp)
	return 0;
    *pitems = (GLshort *) (mp->data);
    *pn = mp->dimensions[0];
    *op = (PyObject *) mp;
    return 1;
}
*/

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

/* not used anywhere
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
*/
#endif				/* Not NUMERIC */

#endif
