#ifdef _PYMOL_WIN32
#include"os_predef.h"
#endif
/**
 *
 * GLU Module for Python (incomplete) 
 * 
 * Modified for incorporation into PyMOL
 *
 * Date: May 99
 *
 * Authors: Thomas Schwaller  <tom.schwaller@linux-magazin.de>
 *          Jim Hugunin       <hugunin@parc.xerox.com>
 *          David Ascher      <da@ski.org>
 *          Michel Sanner     <sanner@scripps.edu>
 * 
 * Notes:
 *
 *   - the file gluconst.py is also needed, as it defines all kinds of
 *     useful constants
 *   - there's a memory leak in the glu tesselator callbacks, but I can't
 *     track it down!
***/

#include "Python.h"

#ifdef CYGWIN
#define WIN32
#endif

#ifdef WIN32
#ifndef MS_WIN32
#define MS_WIN32
#endif
#endif

#ifdef MS_WIN32
#ifndef WIN32
#define WIN32
#endif
#include <windows.h>
#define MCALLBACK (void (__stdcall *)(void))
#define GLUCALLBACK WINAPI
#else
#define MCALLBACK (void (*)(void))
#define GLUCALLBACK
#endif

#ifndef _PYMOL_OSX
#include <GL/gl.h>
#include <GL/glu.h>
#else
#include <gl.h>
#include <glu.h>
#endif


static PyObject *py_gl_Error;

static PyObject *ErrorReturn(char *message)
{
    PyErr_SetString(py_gl_Error, message);
    return NULL;
}
#define ASSERT(E,M) if(!(E)) return ErrorReturn(M)

#define TRY(E) if(! (E)) return NULL

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

typedef struct {
    PyObject_HEAD
    GLUquadricObj * quad;
} gluQuadricObject;

#ifdef GLU_VERSION_1_2
typedef struct {
    PyObject_HEAD
    GLUtesselator * tess;
} gluTesselatorPyObject;

#endif

typedef struct {
    PyObject_HEAD
    PyObject * pyfn;
    void (*fn) (void);
} gluCallbackObject;

/* QUADRIC TYPE STUFF */

static void gluquadric_delete(gluQuadricObject * op)
{
    if (op->quad)
	gluDeleteQuadric(op->quad);
    PyObject_Del(op);
}

static PyObject *gluquadric_getattr(gluQuadricObject * op, char *name)
{
    static struct PyMethodDef py_glu_methods[] =
    {
	{NULL, NULL}
    };
    return Py_FindMethod(py_glu_methods, (PyObject *) op, name);
}

static PyTypeObject GLUquadricType =
{
#ifdef MS_WIN32
    PyObject_HEAD_INIT(NULL)
#else
    PyObject_HEAD_INIT(&PyType_Type)
#endif
    0,				/* Number of items for varobject */
    "glu",			/* Name of this type */
    sizeof(gluQuadricObject),	/* Basic object size */
    0,				/* Item size for varobject */
    (destructor) gluquadric_delete,	/* tp_dealloc */
    0,				/* tp_print */
    (getattrfunc) gluquadric_getattr,	/* tp_getattr */
    0,				/* tp_setattr */
    0,				/* tp_compare */
    0,				/* tp_repr */
    0,				/* tp_as_number */
    0,				/* tp_as_sequence */
    0,				/* tp_as_mapping */
};

#define is_gluQuadricObject(op) ((op)->ob_type == &GLUquadricType)

#define is_gluCallbackObject(op) ((op)->ob_type == &GLUcallbackType)

static GLUquadricObj *getgluquadricvalue(PyObject * op)
#ifdef NUMERIC
{return NULL;};
#else
{
    if (!op || !is_gluQuadricObject(op)) {
	PyErr_BadInternalCall();
	return NULL;
    } else
	return ((gluQuadricObject *) op)->quad;
}
#endif

static PyObject *newgluQuadricObject(GLUquadricObj * Quad)
#ifdef NUMERIC
{return NULL;};
#else
{
    gluQuadricObject *op;
    op = PyObject_NEW(gluQuadricObject, &GLUquadricType);
    if (op == NULL)
	return PyErr_NoMemory();
    op->ob_type = &GLUquadricType;
    op->quad = Quad;
    return (PyObject *) op;
}
#endif

static PyObject *py_glu_NewQuadric(PyObject * self, PyObject * args)
{
    if (!PyArg_ParseTuple(args, ""))
	return NULL;
    return newgluQuadricObject(gluNewQuadric());
}

/* tesselator TYPE STUFF */

#ifdef GLU_VERSION_1_2

static void glutesselator_delete(gluTesselatorPyObject * op)
{
    if (op->tess)
	gluDeleteTess(op->tess);
    PyObject_Del(op);
}

static PyObject *glutesselator_getattr(gluTesselatorPyObject * op, char *name)
{
    static struct PyMethodDef py_glu_methods[] =
    {
	{NULL, NULL}
    };
    return Py_FindMethod(py_glu_methods, (PyObject *) op, name);
}


static PyTypeObject GLUtesselatorType =
{
#ifdef MS_WIN32
    PyObject_HEAD_INIT(NULL)
#else
    PyObject_HEAD_INIT(&PyType_Type)
#endif
    0,				/* Number of items for varobject */
    "gluTesselator",		/* Name of this type */
    sizeof(gluTesselatorPyObject),	/* Basic object size */
    0,				/* Item size for varobject */
    (destructor) glutesselator_delete,	/* tp_dealloc */
    0,				/* tp_print */
    (getattrfunc) glutesselator_getattr,	/* tp_getattr */
    0,				/* tp_setattr */
    0,				/* tp_compare */
    0,				/* tp_repr */
    0,				/* tp_as_number */
    0,				/* tp_as_sequence */
    0,				/* tp_as_mapping */
};

static void suppress_compiler_warnings(void)
{
  gluquadric_delete(NULL);
  gluquadric_getattr(NULL,NULL);
  getgluquadricvalue(NULL);
  newgluQuadricObject(NULL);
  ErrorReturn(NULL);
  suppress_compiler_warnings();
}

#define is_gluTesselatorPyObject(op) ((op)->ob_type == &GLUtesselatorType)
static GLUtesselator *getglutesselatorvalue(PyObject * op)
#ifdef NUMERIC
{return NULL;};
#else
{
    if (!op || !is_gluTesselatorPyObject(op)) {
	PyErr_BadInternalCall();
	return NULL;
    } else
	return ((gluTesselatorPyObject *) op)->tess;
}
#endif

static void GLUCALLBACK beginCB(GLenum type, PyObject * obj)
{
    PyObject *retval = NULL;
    if (obj == NULL)
	return;
    if (!(retval = PyObject_CallMethod(obj, "beginCB", "i", type)))
	PyErr_Print();
    Py_XDECREF(retval);
}

static void GLUCALLBACK endCB(PyObject * obj)
{
    PyObject *retval = NULL;
    if (obj == NULL)
	return;
    if (!(retval = PyObject_CallMethod(obj, "endCB", NULL)))
	PyErr_Print();
    Py_XDECREF(retval);
}

static void GLUCALLBACK errorCB(GLenum errnum, PyObject * obj)
{
    PyObject *retval = NULL;
    if (obj == NULL)
	return;
    if (!(retval = PyObject_CallMethod(obj, "errorCB", "i", errnum)))
	PyErr_Print();
    Py_XDECREF(retval);
}

static void GLUCALLBACK vertexCB(PyObject * vdata, PyObject * obj)
{
    PyObject *retval = NULL;
    if (obj == NULL)
	return;
    Py_INCREF(vdata);
    if (!(retval = PyObject_CallMethod(obj, "vertexCB", "(O)", vdata)))
	PyErr_Print();
    Py_DECREF(vdata);
    Py_XDECREF(retval);
}

static void GLUCALLBACK edgeFlagCB(GLboolean flag, PyObject * obj)
{
    PyObject *retval = NULL;
    if (obj == NULL)
	return;
    if (!(retval = PyObject_CallMethod(obj, "edgeFlagCB", "i", flag)))
	PyErr_Print();
    Py_XDECREF(retval);
}

static void GLUCALLBACK combineCB(GLdouble coords[3],
				  PyObject * vdata,
				  GLfloat weight[4],
				  PyObject ** outData,
				  PyObject * obj)
{
    PyObject *retval;

    if (obj == NULL)
	return;
    retval = PyObject_CallMethod(obj, "combineCB", "((ddd)(ffff)O)",
				 coords[0],
				 coords[1],
				 coords[2],
				 weight[0],
				 weight[1],
				 weight[2],
				 weight[3],
				 vdata);
    if (retval == NULL)
	PyErr_Print();

    Py_XDECREF(retval);
    *outData = retval;
}


static PyObject *py_glu_NewTess(PyObject * self, PyObject * args)
{
    gluTesselatorPyObject *o;
    GLUtesselator *tobj;

    tobj = gluNewTess();
    o = PyObject_NEW(gluTesselatorPyObject, &GLUtesselatorType);
    if (o == NULL)
	return PyErr_NoMemory();
    o->ob_type = &GLUtesselatorType;
    o->tess = tobj;

    /* do not use any of the data-less callbacks */
    gluTessCallback(tobj, GLU_TESS_BEGIN, NULL);
    gluTessCallback(tobj, GLU_TESS_END, NULL);
    gluTessCallback(tobj, GLU_TESS_VERTEX, NULL);
    gluTessCallback(tobj, GLU_TESS_COMBINE, NULL);
    gluTessCallback(tobj, GLU_TESS_EDGE_FLAG, NULL);
    gluTessCallback(tobj, GLU_TESS_ERROR, NULL);

    /* specify our own wrappers for the data-full callbacks */
    gluTessCallback(tobj, GLU_TESS_BEGIN_DATA, MCALLBACK beginCB);
    gluTessCallback(tobj, GLU_TESS_END_DATA, MCALLBACK endCB);
    gluTessCallback(tobj, GLU_TESS_VERTEX_DATA, MCALLBACK vertexCB);
    gluTessCallback(tobj, GLU_TESS_COMBINE_DATA, MCALLBACK combineCB);
    gluTessCallback(tobj, GLU_TESS_EDGE_FLAG_DATA, MCALLBACK(void (*)(void)) edgeFlagCB);
    gluTessCallback(tobj, GLU_TESS_ERROR_DATA, MCALLBACK errorCB);

    return (PyObject *) o;
}



static PyObject *py_glu_TessProperty(PyObject * self, PyObject * args)
{
    PyObject *tess;
    GLdouble value;
    GLenum which;
    GLUtesselator *t;

    if (!PyArg_ParseTuple(args, "Oid", &tess, &which, &value))
	return NULL;
    t = getglutesselatorvalue(tess);
    gluTessProperty(t, which, value);

    Py_INCREF(Py_None);
    return Py_None;
}

#endif

/* cylinder */
static PyObject *py_glu_Cylinder(PyObject * self, PyObject * args)
{
    PyObject *quad;
    GLdouble baseRadius, topRadius, height;
    GLint slices, stacks;
    if (!PyArg_ParseTuple(args, "Odddii", &quad, &baseRadius, &topRadius, &height, &slices, &stacks))
	return NULL;
    gluCylinder(getgluquadricvalue(quad), baseRadius, topRadius, height, slices, stacks);
    Py_INCREF(Py_None);
    return Py_None;
}

/* ---- spheres, partial disks */

static PyObject *py_glu_Sphere(PyObject * self, PyObject * args)
{
    PyObject *quad;
    GLdouble radius;
    GLint slices, stacks;
    if (!PyArg_ParseTuple(args, "Odii", &quad, &radius, &slices, &stacks))
	return NULL;
    gluSphere(getgluquadricvalue(quad), radius, slices, stacks);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_glu_PartialDisk(PyObject * self, PyObject * args)
{
    PyObject *quad;
    GLdouble iradius, oradius, start, sweep;
    GLint slices, loops;
    if (!PyArg_ParseTuple(args, "Oddiidd", &quad, &iradius, &oradius,
			  &slices, &loops, &start, &sweep))
	return NULL;
    gluPartialDisk(getgluquadricvalue(quad), iradius, oradius,
		   slices, loops, start, sweep);
    Py_INCREF(Py_None);
    return Py_None;
}

/* Quadric draw style */

static PyObject *py_glu_QuadricDrawStyle(PyObject * self, PyObject * args)
{
    PyObject *quad;
    GLenum drawStyle;
    if (!PyArg_ParseTuple(args, "Oi", &quad, &drawStyle))
	return NULL;
    gluQuadricDrawStyle(getgluquadricvalue(quad), drawStyle);
    Py_INCREF(Py_None);
    return Py_None;
}

/* Quadric normals */

static PyObject *py_glu_QuadricNormals(PyObject * self, PyObject * args)
{
    PyObject *quad;
    GLenum normals;
    if (!PyArg_ParseTuple(args, "Oi", &quad, &normals))
	return NULL;
    gluQuadricNormals(getgluquadricvalue(quad), normals);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_glu_QuadricOrientation(PyObject * self, PyObject * args)
{
    PyObject *quad;
    GLenum orientation;
    if (!PyArg_ParseTuple(args, "Oi", &quad, &orientation))
	return NULL;
    gluQuadricOrientation(getgluquadricvalue(quad), orientation);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_glu_QuadricTexture(PyObject * self, PyObject * args)
{
    PyObject *quad;
    GLboolean texture;
    if (!PyArg_ParseTuple(args, "Oi", &quad, &texture))
	return NULL;
    gluQuadricTexture(getgluquadricvalue(quad), texture);
    Py_INCREF(Py_None);
    return Py_None;
}


#ifdef GLU_VERSION_1_2

static PyObject *py_glu_TessBeginPolygon(PyObject * self, PyObject * args)
{
    PyObject *tess;
    PyObject *data = NULL;

    if (!PyArg_ParseTuple(args, "O|O", &tess, &data))
	return NULL;

    gluTessBeginPolygon(getglutesselatorvalue(tess), data);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_glu_TessEndPolygon(PyObject * self, PyObject * args)
{
    PyObject *tess;

    if (!PyArg_ParseTuple(args, "O", &tess))
	return NULL;
    gluTessEndPolygon(getglutesselatorvalue(tess));

    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *py_glu_TessBeginContour(PyObject * self, PyObject * args)
{
    PyObject *tess;

    if (!PyArg_ParseTuple(args, "O", &tess))
	return NULL;
    gluTessBeginContour(getglutesselatorvalue(tess));

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_glu_TessEndContour(PyObject * self, PyObject * args)
{
    PyObject *tess;

    if (!PyArg_ParseTuple(args, "O", &tess))
	return NULL;
    gluTessEndContour(getglutesselatorvalue(tess));

    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *py_glu_TessVertex(PyObject * self, PyObject * args)
{
    PyObject *tess;
    PyObject *d = NULL;
    double *vert;
    PyObject *v;
    int n;

    if (!PyArg_ParseTuple(args, "OOO", &tess, &v, &d))
	return NULL;
    TRY(PyArray_AsDoubleArray(&v, &vert, &n));

    Py_INCREF(d);
    gluTessVertex(getglutesselatorvalue(tess), vert, d);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_glu_TessNormal(PyObject * self, PyObject * args)
{
    PyObject *tess;
    double x, y, z;
    if (!PyArg_ParseTuple(args, "Oddd", &tess, &x, &y, &z))
	return NULL;
    gluTessNormal(getglutesselatorvalue(tess), x, y, z);
    Py_INCREF(Py_None);
    return Py_None;
}

#endif				/* GLU_VERSION_1_2 */

static PyObject *py_glu_Perspective(PyObject * self, PyObject * args)
{
    GLdouble a, b, c, d;
    if (!PyArg_ParseTuple(args, "dddd", &a, &b, &c, &d))
	return NULL;
    gluPerspective(a, b, c, d);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_glu_LookAt(PyObject * self, PyObject * args)
{
    GLdouble a, b, c, d, e, f, g, h, i;
    if (!PyArg_ParseTuple(args, "ddddddddd", &a, &b, &c, &d, &e, &f, &g, &h, &i))
	return NULL;
    gluLookAt(a, b, c, d, e, f, g, h, i);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_glu_Ortho2D(PyObject * self, PyObject * args)
{
    GLdouble a, b, c, d;
    if (!PyArg_ParseTuple(args, "dddd", &a, &b, &c, &d))
	return NULL;
    gluOrtho2D(a, b, c, d);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_glu_PickMatrix(PyObject * self, PyObject * args)
{
    int i;
    GLdouble x, y, width, height;
    GLint viewport[4];
    PyObject *viewportlist;
    if (!PyArg_ParseTuple(args, "ddddO!", &x, &y, &width, &height,
			  &PyList_Type, &viewportlist))
	return NULL;
    for (i = 0; i < 4; i++)
	viewport[i] = (int) PyFloat_AsDouble(PyList_GetItem(viewportlist, i));
    gluPickMatrix(x, y, width, height, viewport);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *py_glu_UnProject(PyObject * self, PyObject * args)
{
    GLdouble winx, winy, winz;
    GLdouble objx, objy, objz;
    GLdouble model[16], proj[16];
    GLint view[16];

    if (!PyArg_ParseTuple(args, "ddd", &winx, &winy, &winz))
	return NULL;

    glGetDoublev(GL_MODELVIEW_MATRIX, model);
    glGetDoublev(GL_PROJECTION_MATRIX, proj);
    glGetIntegerv(GL_VIEWPORT, view);
    gluUnProject(winx, winy, winz, model, proj, view, &objx, &objy, &objz);

    return Py_BuildValue("(ddd)", objx, objy, objz);
}

static PyObject *py_glu_Project(PyObject * self, PyObject * args)
{
    GLdouble winx, winy, winz;
    GLdouble objx, objy, objz;
    GLdouble model[16], proj[16];
    GLint view[16];

    if (!PyArg_ParseTuple(args, "ddd", &objx, &objy, &objz))
	return NULL;

    glGetDoublev(GL_MODELVIEW_MATRIX, model);
    glGetDoublev(GL_PROJECTION_MATRIX, proj);
    glGetIntegerv(GL_VIEWPORT, view);
    gluProject(objx, objy, objz, model, proj, view, &winx, &winy, &winz);

    return Py_BuildValue("(ddd)", winx, winy, winz);
}

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


static int typecode2gltype[] =
{-1, GL_UNSIGNED_BYTE, GL_BYTE, GL_SHORT, GL_INT, -1, GL_FLOAT, -1, -1, -1, -1, -1};

static PyObject *
 py_glu_Build2DMipmaps(PyObject * self, PyObject * args)
{
    GLenum target;
    GLint components;
    GLint width;
    GLint height;
    GLenum format;
    GLenum type;
    PyObject *vop;
    PyArrayObject *ap;
    char *data;
    int format_size, type_size, total_size;

    TRY(PyArg_ParseTuple(args, "iiiiiiO", &target, &components, &width, &height, &format, &type, &vop));

    if (PyString_Check(vop)) {
	data = PyString_AsString(vop);
	format_size = glformat_size(format);
	ASSERT(format_size != -1, "invalid format");
	type_size = gltype_size(type);
	ASSERT(format_size != -1, "invalid type");
	total_size = type_size * format_size * width * height / 8;
	ASSERT(total_size >= PyString_Size(vop), "data area too small");

	gluBuild2DMipmaps(target, components, width, height, format, type, data);
    } else {
	if (PyArray_Check(vop)) {
	    ap = (PyArrayObject *) vop;
	    ASSERT(ap->nd == 2 || ap->nd == 3, "array must be either 2 or 3d");
	} else {
	    TRY(ap = (PyArrayObject *) PyArray_ContiguousFromObject(vop, PyArray_UBYTE, 2, 3));
	}

	type = typecode2gltype[ap->descr->type_num];
	ASSERT(type != -1, "can't convert this type of array to an image");

	if (ap->nd == 2) {
	    format = GL_LUMINANCE;
	} else {
	    ASSERT(ap->dimensions[2] == 3 || ap->dimensions[2] == 4, "3d array must be RGB or RGBA");

	    if (ap->dimensions[2] == 3) {
		format = GL_RGB;
	    } else {
		format = GL_RGBA;
	    }
	}
	gluBuild2DMipmaps(target, components, width, height, format, type, ap->data);
    }
    Py_INCREF(Py_None);
    return Py_None;
}
#endif

static PyObject *py_glu_ErrorString(PyObject * self, PyObject * args)
{
    int errCode;

    if (!PyArg_ParseTuple(args, "i", &errCode))
	return NULL;
    return Py_BuildValue("s", gluErrorString((GLenum) errCode));
}

static PyMethodDef py_glu_methods[] =
{
    {"gluNewQuadric", py_glu_NewQuadric, 1},
#ifdef GLU_VERSION_1_2
    {"gluNewTess", py_glu_NewTess, 1},
    {"gluTessProperty", py_glu_TessProperty, 1},
    {"gluTessVertex", py_glu_TessVertex, 1},
    {"gluTessNormal", py_glu_TessNormal, 1},
    {"gluTessBeginPolygon", py_glu_TessBeginPolygon, 1},
    {"gluTessBeginContour", py_glu_TessBeginContour, 1},
    {"gluTessEndContour", py_glu_TessEndContour, 1},
    {"gluTessEndPolygon", py_glu_TessEndPolygon, 1},
#endif
    {"gluSphere", py_glu_Sphere, 1},
    {"gluCylinder", py_glu_Cylinder, 1},
    {"gluPartialDisk", py_glu_PartialDisk, 1},
    {"gluQuadricDrawStyle", py_glu_QuadricDrawStyle, 1},
    {"gluQuadricNormals", py_glu_QuadricNormals, 1},
    {"gluQuadricTexture", py_glu_QuadricTexture, 1},
    {"gluQuadricOrientation", py_glu_QuadricOrientation, 1},
    {"gluPerspective", py_glu_Perspective, 1},
    {"gluLookAt", py_glu_LookAt, 1},
    {"gluOrtho2D", py_glu_Ortho2D, 1},
    {"gluPickMatrix", py_glu_PickMatrix, 1},
    {"gluUnProject", py_glu_UnProject, 1},
    {"gluProject", py_glu_Project, 1},
    {"gluErrorString", py_glu_ErrorString, 1},
#ifdef NUMERIC
    {"gluBuild2DMipmaps", py_glu_Build2DMipmaps, 1},
#endif
    {NULL, NULL}
};


#ifdef NUMERIC

/* C prototype to suppress GCC compiler warning...*/
DL_EXPORT(void)  init_glu_num(void);
DL_EXPORT(void) init_glu_num(void)
#else

/* C prototype to suppress GCC compiler warning...*/
DL_EXPORT(void)  init_glu(void);
DL_EXPORT(void) init_glu(void)
#endif
{
#ifdef MS_WIN32
#ifdef GLU_VERSION_1_2
    GLUtesselatorType.ob_type = &PyType_Type;
#endif
#ifdef NUMERIC
    GLUquadricType.ob_type = &PyType_Type;
#endif
#endif
#ifdef NUMERIC
    Py_InitModule("_glu_num", py_glu_methods);
#else
    Py_InitModule("_glu", py_glu_methods);
#endif
}

/* for distutils compatibility on WIN32 */
#ifndef NUMERIC
/* C prototype to suppress GCC compiler warning...*/
DL_EXPORT(void) init_glumodule(void);
DL_EXPORT(void) init_glumodule(void) {init_glu();}
#endif
