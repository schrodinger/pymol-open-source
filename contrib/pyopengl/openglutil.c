#ifdef _PYMOL_WIN32
#include"os_predef.h"
#endif
/** 
 *
 *  This module provides utility functions to use in conjunction with
 *  the opengl module.  They are packaged separately because they are
 *  not straightforward interfaces to the OpenGL 1.0 or 1.1 library
 *  spec, but "home-grown" functions.
 *
 * WARNING: Modified for usage within PyMOL
 *
**/

#ifdef __MWERKS__
#define NO_IMPORT_ARRAY
#endif

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
#include <GL/glu.h>
#else
#include <gl.h>
#include <glu.h>
#endif

#include <math.h>
#include "Python.h"
#include <string.h>
#include "openglutil.h"

static PyObject *ErrorReturn(char *message)
{
    PyErr_SetString(gl_Error, message);
    return NULL;
}

static void suppress_compiler_warnings(void)
{
  ErrorReturn(NULL);
  suppress_compiler_warnings();
}


/*
 * The procedure gl_SaveTiff was added by 
 * Lothar Birk <lb@ocean.fb12.TU-Berlin.DE>
 * The new method is also listed in the method table 
 * at the end of the file.
 * 
 * Last change: 6. Feb 1997 by lb
 * 
 * gl_SaveTiff saves an PyOpenGL image or a rectangular part of it 
 * into a TIFF file. You must have the TIFF-library (usually 
 * libtiff.a and tiffio.h) at hand to compile and link this procedure.
 * 
 * gl_SaveTiff was derived from gl_save (by Tom Schwaller) and
 * the example code of question "Q_3_18 Saving OpenGL screen output"
 * of the OpenGL-FAQ written by Reto Konradi.
 * 
 * The procedure takes six arguments:
 *      filename      char*       filename for the TIFF-file
 *      Orient        int         orientation flag: 
 *                                  Orient = 0 => save as is
 *                                  Orient = 1 => rotate by 90 degrees
 *      OffX, OffY    int         Pixel coordinates of lower left corner
 *                                  of the rectangle to be saved
 *                                  Use 0,0 to get the complete OpenGL
 *                                  frame.
 *      ImgW          int         width in pixels of image(part) to be saved
 *      ImgH          int         height in pixels of image(part) to be saved
 *
 */

#ifdef LIBTIFF

#include <tiffio.h>


static PyObject *gl_SaveTiff(PyObject * self, PyObject * args)
{
    char *filename;
    int ImgH, ImgW, OffX, OffY, Orient;

    int bufSize, i, res, rowI, rowsPerStrip;
    int tiffH, tiffW;

    GLubyte *buf;
    TIFF *tif;

    TRY(PyArg_ParseTuple(args, "siiiii", &filename, &Orient,
			 &OffX, &OffY,
			 &ImgW, &ImgH));

    tif = TIFFOpen(filename, "w");
    if (tif == NULL) {
	PyErr_SetString(gl_Error, "could not create TIFF file");
	return NULL;
    }
    if (Orient == 0) {
	tiffW = ImgW;
	tiffH = ImgH;
	bufSize = 4 * ((3 * tiffW + 3) / 4);
	glPixelStorei(GL_PACK_ALIGNMENT, 4);
    } else {
	tiffW = ImgH;
	tiffH = ImgW;
	bufSize = 3 * tiffW;
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
    }

    rowsPerStrip = (8 * 1024) / (3 * tiffW);
    if (rowsPerStrip == 0)
	rowsPerStrip = 1;

    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, tiffW);
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH, tiffH);
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_LZW);
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
    TIFFSetField(tif, TIFFTAG_FILLORDER, FILLORDER_MSB2LSB);
    TIFFSetField(tif, TIFFTAG_DOCUMENTNAME, filename);
    TIFFSetField(tif, TIFFTAG_IMAGEDESCRIPTION, "PyOpenGL screen dump");
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 3);
    TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, rowsPerStrip);
    TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);

    buf = PyMem_NEW(GLubyte, bufSize);

    res = 0;
    for (rowI = 0; rowI < tiffH; rowI++) {
	if (Orient == 0)
	    glReadPixels(OffX, ImgH - 1 - rowI + OffY, ImgW, 1,
			 GL_RGB, GL_UNSIGNED_BYTE, buf);
	else
	    glReadPixels(rowI + OffX, OffY, 1, ImgH,
			 GL_RGB, GL_UNSIGNED_BYTE, buf);

	if (TIFFWriteScanline(tif, buf, rowI, 0) < 0) {
	    PyErr_SetString(gl_Error, "error while writing TIFF file");
	    res = 1;
	    break;
	}
    }

    PyMem_DEL(buf);

    TIFFFlushData(tif);
    TIFFClose(tif);

    if (res == 0) {
	Py_INCREF(Py_None);
	return Py_None;
    } else
	return NULL;
}

#endif				/* LIBTIFF */

static PyObject *gl_SavePPM(PyObject * self, PyObject * args)
{
    char *name;
    int i, width, height;
    GLubyte *pixelbuffer;
    FILE *fp;

    if (!PyArg_ParseTuple(args, "sii", &name, &width, &height))
	return NULL;

    pixelbuffer = PyMem_NEW(GLubyte, 3 * width * height);
    glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixelbuffer);
    fp = fopen(name, "wb");

    if (!fp) {
	PyErr_SetString(PyExc_IOError, "error while opening file");
	return NULL;
    }
    fprintf(fp, "P6\n# Python OpenGL\n%d %d\n255\n", width, height);

    for (i = height - 1; i >= 0; i--)
	fwrite(pixelbuffer + i * width * 3, 1, width * 3, fp);

    fclose(fp);
    PyMem_DEL(pixelbuffer);

    Py_INCREF(Py_None);
    return Py_None;
}

/*
 * Generate EPS file.
 * Contributed by Miguel A. De Riera Pasenau <miguel@DALILA.UPC.ES>
 */

static PyObject *gl_SaveEPS(PyObject * self, PyObject * args)
{
    char *name;
    int i, width, height, pos, components = 3;
    GLubyte *pixels;
    FILE *fp;
    unsigned char *curpix, bitpixel;
    double pix;

    if (!PyArg_ParseTuple(args, "sii", &name, &width, &height))
	return NULL;

    pixels = PyMem_NEW(GLubyte, 3 * width * height);

    glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixels);

    fp = fopen(name, "wb");

    if (!fp) {
	PyErr_SetString(PyExc_IOError, "error while opening file");
	return NULL;
    }
    fprintf(fp, "%%!PS-Adobe-2.0 EPSF-1.2\n");
    fprintf(fp, "%%%%Creator: OpenGL pixmap render output\n");
    fprintf(fp, "%%%%BoundingBox: 0 0 %d %d\n", width, height);
    fprintf(fp, "%%%%EndComments\n");

    i = (((width * height) + 7) / 8) / 40;	/* # of lines, 40 bytes per line */
    fprintf(fp, "%%%%BeginPreview: %d %d %d %d\n%%", width, height, 1, i);
    pos = 0;
    curpix = (unsigned char *) pixels;
    for (i = 0; i < width * height * components;) {
	bitpixel = 0;
	pix = 0.0;
	pix = 0.30 * (double) curpix[i++] + 0.59;
   pix += (double) curpix[i++];
   pix += 0.11 * (double) curpix[i++];
	if (pix > 127.0)
	    bitpixel |= 0x80;
	pix = 0.30 * (double) curpix[i++] + 0.59;
   pix += (double) curpix[i++];
   pix += 0.11 * (double) curpix[i++];
	if (pix > 127.0)
	    bitpixel |= 0x40;
	pix = 0.30 * (double) curpix[i++] + 0.59;
   pix += (double) curpix[i++];
   pix += 0.11 * (double) curpix[i++];
	if (pix > 127.0)
	    bitpixel |= 0x20;
	pix = 0.30 * (double) curpix[i++] + 0.59;
   pix += (double) curpix[i++];
   pix += 0.11 * (double) curpix[i++];
	if (pix > 127.0)
	    bitpixel |= 0x10;
	pix = 0.30 * (double) curpix[i++] + 0.59;
   pix += (double) curpix[i++];
   pix += 0.11 * (double) curpix[i++];
	if (pix > 127.0)
	    bitpixel |= 0x08;
	pix = 0.30 * (double) curpix[i++] + 0.59;
   pix += (double) curpix[i++];
   pix += 0.11 * (double) curpix[i++];
	if (pix > 127.0)
	    bitpixel |= 0x04;
	pix = 0.30 * (double) curpix[i++] + 0.59;
   pix += (double) curpix[i++];
   pix += 0.11 * (double) curpix[i++];
	if (pix > 127.0)
	    bitpixel |= 0x02;
	pix = 0.30 * (double) curpix[i++] + 0.59;
   pix += (double) curpix[i++];
   pix += 0.11 * (double) curpix[i++];
	if (pix > 127.0)
	    bitpixel |= 0x01;

	fprintf(fp, "%02hx", bitpixel);
	if (++pos >= 40) {
	    fprintf(fp, "\n%%");
	    pos = 0;
	}
    }
    if (pos)
	fprintf(fp, "\n%%%%EndPreview\n");
    else
	fprintf(fp, "%%EndPreview\n");

    fprintf(fp, "gsave\n");
    fprintf(fp, "/bwproc {\n");
    fprintf(fp, "    rgbproc\n");
    fprintf(fp, "    dup length 3 idiv string 0 3 0\n");
    fprintf(fp, "    5 -1 roll {\n");
    fprintf(fp, "    add 2 1 roll 1 sub dup 0 eq\n");
    fprintf(fp, "    { pop 3 idiv 3 -1 roll dup 4 -1 roll dup\n");
    fprintf(fp, "        3 1 roll 5 -1 roll put 1 add 3 0 }\n");
    fprintf(fp, "    { 2 1 roll } ifelse\n");
    fprintf(fp, "    } forall\n");
    fprintf(fp, "    pop pop pop\n");
    fprintf(fp, "} def\n");
    fprintf(fp, "systemdict /colorimage known not {\n");
    fprintf(fp, "    /colorimage {\n");
    fprintf(fp, "        pop\n");
    fprintf(fp, "        pop\n");
    fprintf(fp, "        /rgbproc exch def\n");
    fprintf(fp, "        { bwproc } image\n");
    fprintf(fp, "    } def\n");
    fprintf(fp, "} if\n");
    fprintf(fp, "/picstr %d string def\n", width * components);
    fprintf(fp, "%d %d scale\n", width, height);
    fprintf(fp, "%d %d %d\n", width, height, 8);
    fprintf(fp, "[%d 0 0 %d 0 0]\n", width, height);
    fprintf(fp, "{currentfile picstr readhexstring pop}\n");
    fprintf(fp, "false %d\n", components);
    fprintf(fp, "colorimage\n");

    curpix = (unsigned char *) pixels;
    pos = 0;
    for (i = width * height * components; i > 0; i--) {
	fprintf(fp, "%02hx", *curpix++);
	if (++pos >= 40) {
	    fprintf(fp, "\n");
	    pos = 0;
	}
    }
    if (pos)
	fprintf(fp, "\n");

    fprintf(fp, "grestore\n");

    fclose(fp);
    PyMem_DEL(pixels);

    Py_INCREF(Py_None);
    return Py_None;
}

/* 
   Picking function 
   If you want the more traditional 0.0 to 1.0 numbers, convert to a long, 
   add sys.maxint, then divide by (sys.maxint*2).
   If you want the physical depth, multiply that by the frustrum depth and 
   add your near clipping plane.

   Released to the public domain, Mike Fletcher 1999
 */

static PyObject *processHits(GLint hits, GLuint buffer[])
{
    /* convert a select buffer into a set of tuples
       for easy processing on the Python side
       ((minz, maxz, (names,...)), (minz, maxz, (names,...)),... ) */
    PyObject *base;
    PyObject *current;
    PyObject *nameset;
    int i, j, names;
    GLuint *ptr;

    base = PyTuple_New((int) hits);
    ptr = (GLuint *) buffer;

    for (i = 0; i < hits; i++) {	/* for each hit  */
	names = (int) *ptr;	/* names is now dereferenced to number of names */
	ptr++;
	/* build the current tuple */
	current = PyTuple_New(3);	/* min, max, names */
	/* add current to the base set */
	PyTuple_SetItem(base, i, current);
	/* fill in min and max */
	PyTuple_SetItem(current, 0, PyLong_FromUnsignedLong(*ptr));
	ptr++;
	PyTuple_SetItem(current, 1, PyLong_FromUnsignedLong(*ptr));
	ptr++;
	/* add each name to a new tuple that's the third element */
	nameset = PyTuple_New(names);
	PyTuple_SetItem(current, 2, nameset);
	for (j = 0; j < names; j++) {	/* for each name */
	    PyTuple_SetItem(nameset, j, PyInt_FromLong(*ptr));
	    ptr++;
	}
    }
    return base;
}
static char gl_SelectWithCallback__doc__[] =
"glSelectWithCallback(int x, int y, Callable callback, int xsize=5, int ysize=5)\n\
  x,y -- x and y window coordinates for the center of the pick box\n\
  rendercallback -- callback (callable Python object) taking 0 arguments\n\
    which performs pick-mode rendering\n\
  xsize,ysize -- x and y dimensions of the pick box (default = 5x5)\n\
\n\
The function returns a tuple (possibly empty) of:\n\
  ( (minimumzdepth, maximumzdepth, (name, name, name),...)\n\
    minimumzdepth, maximumzdepth -- depths in integer format\n\
      If you want the more traditional 0.0 to 1.0 numbers, divide\n\
	  by (2**32)-1\n\
      If you want the physical depth, multiply that by the frustrum depth and\n\
        add your near clipping plane.\n\
    name -- the names (integers) used in calls to glPushName( int )";

static PyObject *gl_SelectWithCallback(PyObject * self, PyObject * args)
{
    /* perform traditional OpenGL selection (i.e. using render mode GL_SELECT
       code is adapted from the red book).  Note: you should set up your viewing
       matrix _before_ calling this function, as it just uses the current projection
       matrix.  Callback should take 0 arguments.

       glSelectWithCallback( x, y, rendercallback, xsize, ysize );
       x,y -- x and y window coordinates for the center of the pick box
       rendercallback -- callback taking 0 arguments which performs the rendering
       xsize,ysize -- x and y dimensions of the pick box (default = 5x5)
     */
    int x, y;			/* the pick coordinates */
    PyObject *rendercallback;	/* object to call to perform pick-rendering */
    double xsize = 5;
    double ysize = 5;		/* horizontal and vertical size of pick window */

    GLuint selectBuffer[512];	/* up to 512 hits per call by default, no idea what happens if this goes kablooie, hopefully OpenGL does internal checking... */
    GLint hits;
    GLint viewport[4];
    GLdouble previousviewmatrix[16];

    PyObject *result;		/* result to return */

    if (!(PyArg_ParseTuple(args, "iiO|dd", &x, &y, &rendercallback, &xsize, &ysize))) {
	return NULL;
    }
    /* don't bother if the callback isn't going to work... */
    if (PyCallable_Check(rendercallback) == 1) {
	/* query to determine the current viewport size */
	glGetIntegerv(GL_VIEWPORT, viewport);

	/* set selectBuffer as our current selection buffer */
	glSelectBuffer(512, selectBuffer);
	/* switch to selection rendering */
	(void) glRenderMode(GL_SELECT);

	/* initialise name stack */
	glInitNames();
	/* setup the pick matrix */
	glMatrixMode(GL_PROJECTION);
	glGetDoublev(GL_PROJECTION_MATRIX, previousviewmatrix);
	/*glPushMatrix(); */
	glLoadIdentity();
	gluPickMatrix((GLdouble) x, (GLdouble) (viewport[3] - y),
		      (GLdouble) xsize, (GLdouble) ysize, viewport);
	glMultMatrixd(previousviewmatrix);	/* we want to transform this by pick, not other way around? */
	/* we already decided we could call this object... */
	PyObject_CallObject(rendercallback, NULL);
	/* flush the GL, so that we know everything has been rendered */
	glFlush();
	glMatrixMode(GL_PROJECTION);
	/*glPopMatrix(); */
	glLoadMatrixd(previousviewmatrix);

	/* switch mode, getting a count of hits */
	hits = glRenderMode(GL_RENDER);
	/* do the processing */
	result = processHits(hits, selectBuffer);
    } else {
	result = PyTuple_New(0);
    }
    Py_INCREF(result);
    return result;
}

/* experimental picking function ends */
/* experimental */

static PyObject *gl_ColorVertex2d(PyObject * self, PyObject * args)
{
    GLdouble arg1, arg2, *vert;
    PyObject *op1;
    int n;
#ifdef NUMERIC
    GLdouble *x, *y, *c;
    PyObject *op2, *colors = NULL;
    PyArrayObject *cp = NULL;
    int m, i, ok1, ok2, nc;
#endif
    if (PyArg_ParseTuple(args, "dd", &arg1, &arg2))
	glVertex2d(arg1, arg2);
    else {
	PyErr_Clear();
	if (PyArg_ParseTuple(args, "O", &op1)) {
	    TRY(PyArray_AsDoubleArray(&op1, &vert, &n));
	    if (n < 2) {
		PyErr_SetString(gl_Error, "need element with at least 2 items");
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
		PyErr_SetString(gl_Error, "coordinate arrays must be of same length or not enough memory");
		PyArray_ClearMemory(op1, x);
		PyArray_ClearMemory(op2, y);
		return NULL;
	    }
	    if (colors) {
		if (!(cp = (PyArrayObject *) PyArray_ContiguousFromObject(colors,
						PyArray_DOUBLE, 1, 2))) {
		    PyArray_ClearMemory(op1, x);
		    PyArray_ClearMemory(op2, y);
		    return NULL;
		}
		c = (GLdouble *) cp->data;
		nc = PyArray_Size((PyObject *) cp);
		if (((nc % 3) != 0) || (n != nc / 3)) {
		    PyErr_SetString(gl_Error, "wrong color matrix size");
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

#ifdef NUMERIC
/* MS March 4 '99 */

/********************************************************************
  Tries to create a contiguous numeric array of type typecode from a 
  Python object. Works for list, tuples and numeric arrays.

  obj: Numeric array Python object
  typecode: data type PyArray_{ CHAR, UBYTE, SBYTE, SHORT, INT, LONG, FLOAT,
                                DOUBLE, CFLOAT, CDOUBLE }
  expectnd: required number of dimensions. Used for checking. Ignored if <=0.
  expectdims: array of expected extends. Used for checking. Ignored if <=0.

  Raises ValueError exceptions if:
  - the PyArray_ContiguousFromObject fails
  - the array has a bad shape
  - the extent of a given dimension doesn't match the specified extent.
********************************************************************/

static PyArrayObject *contiguous_typed_array(PyObject * obj, int typecode,
					   int expectnd, int *expectdims)
{
    PyArrayObject *arr;
    int i;
    char buf[255];

    /* if the shape and type are OK, this function increments the reference
       count and arr points to obj */
    if ((arr = (PyArrayObject *) PyArray_ContiguousFromObject(obj,
							      typecode, 0,
							  10)) == NULL) {
	sprintf(buf, "Failed to make a contiguous array of type %d\n", typecode);
	PyErr_SetString(PyExc_ValueError, buf);
	return NULL;
    }
    if (expectnd > 0) {
	if (arr->nd > expectnd + 1 || arr->nd < expectnd) {
	    Py_DECREF((PyObject *) arr);
	    PyErr_SetString(PyExc_ValueError,
			    "Array has wrong number of dimensions");
	    return NULL;
	}
	if (arr->nd == expectnd + 1) {
	    if (arr->dimensions[arr->nd - 1] != 1) {
		Py_DECREF((PyObject *) arr);
		PyErr_SetString(PyExc_ValueError,
				"Array has wrong number of dimensions");
		return NULL;
	    }
	}
	if (expectdims) {
	    for (i = 0; i < expectnd; i++)
		if (expectdims[i] > 0)
		    if (expectdims[i] != arr->dimensions[i]) {
			Py_DECREF((PyObject *) arr);
			sprintf(buf, "The extent of dimension %d is %d while %d was expected\n",
				i, arr->dimensions[i], expectdims[i]);
			PyErr_SetString(PyExc_ValueError, buf);
			return NULL;
		    }
	}
    }
    return arr;
}

/*
   When composing successively rotations the resulting matrix often
   becomes nonorthogonal leading to skewing and scaling effects.
   This function takes a 4x4 matrix of DOUBLE that represents a OpenGL
   transformation and reorthogonalizes it
   uses: contiguous_typed_array
   available from the interpreter as:
   CleanRotMat( mat )
 */
static PyObject *gl_CleanRotMat(PyObject * self, PyObject * args)
{
    PyArrayObject *outmat, *inmat;
    PyObject *matrix;
    double a[4][4];
    int expected_dims[2] =
    {4, 4};

    float s;
    int i;

    if (!PyArg_ParseTuple(args, "O", &matrix))
	return NULL;

    inmat = contiguous_typed_array(matrix, PyArray_DOUBLE, 2, expected_dims);
    if (!inmat)
	return NULL;

    memcpy(a, (void *) inmat->data, sizeof(double) * 16);

    for (i = 0; i < 3; i++)
	a[i][3] = a[3][i] = 0.0;
    a[3][3] = 1.0;

    for (i = 0, s = 0.0; i < 3; i++)
	s += (float)(a[0][i] * a[0][i]);
    s = (float)sqrt(s);
    for (i = 0; i < 3; i++)
	a[0][i] /= s;		/* first row normalized */

    a[2][0] = a[0][1] * a[1][2] - a[0][2] * a[1][1];
    a[2][1] = a[0][2] * a[1][0] - a[0][0] * a[1][2];
    a[2][2] = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    for (i = 0, s = 0.0; i < 3; i++)
	s += (float)(a[2][i] * a[2][i]);
    s = (float)sqrt(s);
    for (i = 0; i < 3; i++)
	a[2][i] /= s;		/* third row orthonormal to first */

    a[1][0] = a[2][1] * a[0][2] - a[2][2] * a[0][1];
    a[1][1] = a[2][2] * a[0][0] - a[2][0] * a[0][2];
    a[1][2] = a[2][0] * a[0][1] - a[2][1] * a[0][0];
    for (i = 0, s = 0.0; i < 3; i++)
	s += (float)(a[1][i] * a[1][i]);
    s = (float)sqrt(s);
    for (i = 0; i < 3; i++)
	a[1][i] /= s;		/* second row orthonormal to 1,3 */

    outmat = (PyArrayObject *) PyArray_FromDims(2, expected_dims, PyArray_DOUBLE);
    if (!outmat) {
	PyErr_SetString(PyExc_RuntimeError,
			"Failed to allocate memory for matrix");
	return NULL;
    }
    memcpy(outmat->data, a, 16 * sizeof(double));
    return (PyObject *) outmat;
}

/*************
  Compute vector v, normal to the triangle (p1,p2,p3) assuming that the order
  of the points p1,p2,p3 provides the face's orientation
**************/
static void triangle_normal(double *p1, double *p2, double *p3, float *v)
{
    double v1[3], v2[3], norm;
    short i;

    for (i = 0; i < 3; i++) {
	v1[i] = p2[i] - p1[i];	/* vector (p1,p2) */
	v2[i] = p3[i] - p2[i];	/* vector (p2,p3) */
    }
    v[0] = (float)(v1[1] * v2[2] - v1[2] * v2[1]);	/* v3 = v1^v2 */
    v[1] = (float)(v1[2] * v2[0] - v1[0] * v2[2]);
    v[2] = (float)(v1[0] * v2[1] - v1[1] * v2[0]);

    norm = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (norm > 0.01) {
	for (i = 0; i < 3; i++)
	    v[i] = (float)(v[i]/norm);
    } else {
	for (i = 0; i < 3; i++)
	    v[i] = 0.0;
    }
}

#include <string.h>
/*
   Compute normal vector for a bunch of triangles specified as:
   an n*3 sequence of floats for the vertices coordinates and
   an m*3 array of integers for the triangle's topology
   The thirds and optional (string) argument specifies the computation mode.
   This function can work in three different modes (default: 'PER_VERTEX'):
   'PER_FACE': computes the vector normal to each triangle using
   triangle_normal. The resulting normals are returned
   in a m*3 array of floats.
   'PER_VERTEX': after the face normals have been computed they normals
   for each vertex are obtained by summing up the faces
   normals of each triangle this vertex belongs to. The 
   resulting normals are returned in a n*3 array of floats.
   'BOTH': The face and vertex normals are computed and both are returned.

   uses: contiguous_typed_array, triangle_normal
   available from the interpreter as:
   TriangleNormals( vertices, triangles, | mode )
 */

static PyObject *gl_TriangleNormals(PyObject * self, PyObject * args)
{
    PyObject *vert, *tri;
    double *v_data;
    int *t_data;
    int expected_dims[2] =
    {0, 3};			/* expected dimensions */
    PyArrayObject *out, *out1, *inv, *intr;
    char *s = "PER_FACE";
    int i, j, k, *tric;
    float *trinorm, *vnorm;
    char buf[255];

    TRY(PyArg_ParseTuple(args, "OO|s", &vert, &tri, &s));

    /* check the shape of the vertex sequence and get a handle to the data */
    inv = contiguous_typed_array(vert, PyArray_DOUBLE, 2, expected_dims);
    if (!inv)
	return NULL;
    v_data = (double *) inv->data;

    /* check the shape of the triangles sequence and get a handle to the data */
    intr = contiguous_typed_array(tri, PyArray_INT, 2, expected_dims);
    if (!intr)
	return NULL;
    t_data = (int *) intr->data;

    /* compute the triangle's normals */
    trinorm = (float *) malloc(intr->dimensions[0] * 3 * sizeof(float));
    if (!trinorm) {
	PyErr_SetString(PyExc_RuntimeError,
			"Failed to allocate memory for the normals");
	return NULL;
    }
    for (i = 0; i < 3 * intr->dimensions[0]; i += 3) {
	int v1, v2, v3;
	v1 = t_data[i];
	if (v1 >= inv->dimensions[0]) {
	    sprintf(buf, "Coordinates index %d in face %d out of range", v1,
		    inv->dimensions[0]);
	    PyErr_SetString(PyExc_ValueError, buf);
	    return NULL;
	}
	v2 = t_data[i + 1];
	if (v2 >= inv->dimensions[0]) {
	    sprintf(buf, "Coordinates index %d in face %d out of range", v2,
		    inv->dimensions[0]);
	    PyErr_SetString(PyExc_ValueError, buf);
	    return NULL;
	}
	v3 = t_data[i + 2];
	if (v3 >= inv->dimensions[0]) {
	    sprintf(buf, "Coordinates index %d in face %d out of range", v3,
		    inv->dimensions[0]);
	    PyErr_SetString(PyExc_ValueError, buf);
	    return NULL;
	}
	triangle_normal(&v_data[v1 * 3], &v_data[v2 * 3], &v_data[v3 * 3],
			&trinorm[i]);
    }

    if (strncmp(s, "PER_FACE", strlen(s)) != 0)
	/* compute the vertices normals */
    {
	vnorm = (float *) malloc(inv->dimensions[0] * 3 * sizeof(float));
	tric = (int *) malloc(inv->dimensions[0] * sizeof(int));
	if (!vnorm || !tric) {
	    PyErr_SetString(PyExc_RuntimeError,
			    "Failed to allocate memory for the normals");
	    return NULL;
	}
	for (i = 0; i < inv->dimensions[0]; i++) {
	    tric[i] = 0;
	    for (j = 0; j < 3; j++)
		vnorm[i * 3 + j] = 0.0;
	}
	for (i = 0; i < intr->dimensions[0] * 3; i += 3) {	/* loop over triangles */
	    for (k = 0; k < 3; k++) {	/* loop over vertices */
		tric[t_data[i + k]]++;
		vnorm[3 * t_data[i + k]] += trinorm[i];
		vnorm[3 * t_data[i + k] + 1] += trinorm[i + 1];
		vnorm[3 * t_data[i + k] + 2] += trinorm[i + 2];
	    }
	}
	for (i = 0; i < inv->dimensions[0]; i++) {
	    for (k = 0; k < 3; k++)
		vnorm[i * 3 + k] /= tric[i];
	}
	free(tric);

	out = (PyArrayObject *) PyArray_FromDimsAndData(2, inv->dimensions,
					  PyArray_FLOAT, (char *) vnorm);
	if (!out) {
	    PyErr_SetString(PyExc_RuntimeError,
			    "Failed to allocate memory for normals");
	    return NULL;
	}
	out->flags |= OWN_DATA;	/* so we'll free this memory when this 
				   array will be garbage collected */

	if (strncmp(s, "BOTH", strlen(s)) == 0) {
	    out1 = (PyArrayObject *) PyArray_FromDimsAndData(2, intr->dimensions,
					PyArray_FLOAT, (char *) trinorm);
	    if (!out1) {
		PyErr_SetString(PyExc_RuntimeError,
				"Failed to allocate memory for normals");
		return NULL;
	    }
	    out1->flags |= OWN_DATA;	/* so we'll free this memory when this 
					   array will be garbage collected */

	    return Py_BuildValue("OO", (PyObject *) out, (PyObject *) out1);

	} else {
	    free(trinorm);
	    return (PyObject *) out;
	}
    } else {
	out = (PyArrayObject *) PyArray_FromDimsAndData(2, intr->dimensions,
					PyArray_FLOAT, (char *) trinorm);
	if (!out) {
	    PyErr_SetString(PyExc_RuntimeError,
			    "Failed to allocate memory for normals");
	    return NULL;
	}
	out->flags |= OWN_DATA;	/* so we'll free this memory when this 
				   array will be garbage collected */
	return (PyObject *) out;
    }
}

/*
   WARNING this function is still experimental and has not yet being fully
   tested!
   Build a displaylist for a set of indexed GL geometries (all but GL_POINTS).
   Indexed geometries are specified using a sequence of vertices (n*3 floats)
   and a sequence on m*p integers specifying how to connect the vertices.
   - All parts (faces or lines) can have exactly length p (p vertices) or can
   be of length < p in which case the list is terminated by -1.
   - Sets of lines of length 2, sets of triangles and sets of quads are drawn
   without sending glBegin( ... ) / glEnd () for each primitive. This is more
   efficient but has the draw back that the picking mechanism can no more
   distinguish between different parts. For instance, in order to get 
   pickable triangles one should use GL_POLYGON as the first argument rather
   than GL_TRIANGLES.
   - All arguments are named and some are optional.
   - Lighting is enable only if normals are specified and their number is
   either 1, number of index lists or number of vertices. The normal binding
   mode (OVERALL, PER_PART, PER_VERTEX) is infered from the number of 
   normals available.
   - When no normals are given, lighting is disabled and the polygon mode
   is set to lines.
   - Front and back polygon properties can be specified optionally. They have
   to be a sequence of 5 sequences of properties: one each of the material
   property component: GL_AMBIENT, GL_DIFFUSE, GL_SPECULAR, GL_EMISSION,
   GL_SHININESS
   * The first four property lists have to be (possibly empty) 4-sequences of
   floating point values ranging from 0.0 to 1.0 (RGBA).
   * The fifth property (GL_SHININESS) is a (possibly empty) sequence of
   floating points values ranging from 0.0 to 128.0.
   * property binding modes are selected automatically and separately for
   each property, based on the number of properties available and the number
   of vertices and parts (see normals bind mode).

   Required arguments:
   type       :    GL_LINES, GL_LINE_STRIP, GL_LINE_LOOP, GL_POLYGON,
   GL_QUADS, GL_QUAD_STRIP, GL_TRIANGLES, GL_TRIANGLE_STRIP,
   GL_TRIANGLE_FAN
   coordinates:    sequence of shape n*3 for n vertices (float or double)
   indices    :    sequence of shape m*3 for m triangles (int)

   Optional arguments:
   normals:        sequence of p*3 (float or double).
   If p==n normals will be bound by vertex
   If p==m normals will be bound by face
   If p==1 normal set OVERALL
   else lighting is turned off, GL_LINE is used ?
   frontMaterial:  sequence of p*4 float (0.0-1.0)
   backMaterial:   sequence of p*4 float (0.0-1.0)
   frontMatBind:   binding mode for front material
   backMatBind:    binding mode for back material
   frontAndBack:   int = 0 or 1 to use the front properties for back facing
   polygons.
   texIndices:     texture indices n * 1,2,3 or 4 for n vertices 

   03 '99: Remove overall coloring from display list (has to be done outside)
   03 '99: Added color memory to minimize context switches

   TODO: materials should be indexed too
 */
static short isNewColor(float *c)
{
    static float col[4];
    if (!c) {
	col[0] = col[1] = col[2] = col[3] = -1.0;
	return 0;
    } else {
	if (fabs(c[0] - col[0]) < 0.0001 && fabs(c[1] - col[1]) < 0.0001 &&
	  fabs(c[2] - col[2]) < 0.0001 && fabs(c[3] - col[3]) < 0.0001) {
	    return 0;
	}
#ifdef DEBUG
	printf("new color %f %f %f %f\n", c[0], c[1], c[2], c[3]);
#endif
	col[0] = c[0];
	col[1] = c[1];
	col[2] = c[2];
	col[3] = c[3];
	return 1;
    }
}

static short isNewMaterial(int face, int prop, float *c)
{
    static float col[2][5][4];
    int f, i, j, k;

    if (!c) {
	for (i = 0; i < 2; i++)
	    for (j = 0; j < 5; j++)
		for (k = 0; k < 4; k++)
		    col[i][j][k] = -1.0;
	return 0;
    } else {
	f = (face == GL_FRONT) ? 0 : 1;
	if (fabs(c[0] - col[f][prop][0]) < 0.0001 &&
	    fabs(c[1] - col[f][prop][1]) < 0.0001 &&
	    fabs(c[2] - col[f][prop][2]) < 0.0001 &&
	    fabs(c[3] - col[f][prop][3]) < 0.0001) {
	    return 0;
	}
#ifdef DEBUG
	printf("new material %d %d %f %f %f %f\n", face, prop, c[0], c[1], c[2], c[3]);
#endif
	col[f][prop][0] = c[0];
	col[f][prop][1] = c[1];
	col[f][prop][2] = c[2];
	col[f][prop][3] = c[3];
	return 1;
    }
}

static PyObject *gl_indexedGeomDSPL(PyObject * self, PyObject * args,
				    PyObject * kw)
{
    static char *argnames[] =
    {"type", "coordinates", "indices", "normals",
     "texIndices", "frontMaterial", "backMaterial",
     "frontMatBind", "backMatBind",
     "frontAndBack", "noLightCol",
     NULL};

    int type, frontAndBack = 0, noCol = 0;
    PyObject *py_coords, *py_texInd = NULL, *py_indices, *py_norm = NULL,
    *py_frontMat = NULL, *py_backMat = NULL, *py_frontMatBind = NULL,
    *py_backMatBind = NULL, *lstitem;
    PyArrayObject *c_array, *ind_array, *tex_array = NULL, *n_array = NULL,
    *fmb_array = NULL, *bmb_array = NULL, *fm_array[5], *bm_array[5];

    double *coords, *norm = NULL;
    float *frontMat[5], *backMat[5], *texInd = NULL;

    int expected_coord_dims[2] =
    {-1, 3}, expected_prop_dims[2] =
    {-1, 4};

    int *indices;
    int expected_bind_dims = 5, hasFrontMatBind = 0, hasBackMatBind = 0;

    GLuint dpl;
    int i, j, k, l, v, fixed = 0, face, normBinding, frontMatBind[5],
     backMatBind[5], *fmb, *bmb, NONE = -1, OVERALL = 10, PER_VERTEX = 11,
     PER_PART = 12, propConst[] =
    {GL_AMBIENT, GL_DIFFUSE, GL_SPECULAR, GL_EMISSION,
     GL_SHININESS};
    char buf[255];

    if (!PyArg_ParseTupleAndKeywords(args, kw, "iOO|OOOOOOii", argnames,
				     &type, &py_coords, &py_indices,
				     &py_norm, &py_texInd,
				     &py_frontMat, &py_backMat,
				     &py_frontMatBind, &py_backMatBind,
				     &frontAndBack, &noCol))
	return NULL;

    c_array = contiguous_typed_array(py_coords, PyArray_DOUBLE, 2,
				     expected_coord_dims);
    if (!c_array)
	return NULL;
    coords = (double *) c_array->data;

    ind_array = contiguous_typed_array(py_indices, PyArray_INT, 2, NULL);
    if (!ind_array) {
	Py_DECREF((PyObject *) c_array);
	return NULL;
    }
    indices = (int *) ind_array->data;

    if (py_texInd == Py_None) {
	Py_DECREF(Py_None);
	py_texInd = NULL;
    }
    if (py_texInd) {
	tex_array = contiguous_typed_array(py_texInd, PyArray_FLOAT, 2, NULL);
	if (!tex_array) {
	    Py_DECREF((PyObject *) c_array);
	    Py_DECREF((PyObject *) ind_array);
	    return NULL;
	}
	texInd = (float *) tex_array->data;

	if (tex_array->dimensions[0] != c_array->dimensions[0]) {
	    PyErr_SetString(PyExc_RuntimeError,
	    "Number of texture indices doesn't match number of vertices");
	    Py_DECREF((PyObject *) c_array);
	    Py_DECREF((PyObject *) ind_array);
	    Py_DECREF((PyObject *) tex_array);
	    return NULL;
	}
    }
    if (py_norm == Py_None) {
	Py_DECREF(Py_None);
	py_norm = NULL;
    }
    if (py_norm) {
	n_array = contiguous_typed_array(py_norm, PyArray_DOUBLE, 2, NULL);
	if (!n_array) {
	    Py_DECREF((PyObject *) c_array);
	    Py_DECREF((PyObject *) ind_array);
	    if (tex_array) {
		Py_DECREF((PyObject *) tex_array);
	    }
	    return NULL;
	}
	norm = (double *) n_array->data;

	if (n_array->dimensions[0] == c_array->dimensions[0])
	    normBinding = PER_VERTEX;
	else if (n_array->dimensions[0] == ind_array->dimensions[0])
	    normBinding = PER_PART;
	else if (n_array->dimensions[0] == 1)
	    normBinding = OVERALL;
	else
	    normBinding = NONE;
    } else
	normBinding = NONE;

    /* check front material binding parameter */
    if (py_frontMatBind == Py_None) {
	Py_DECREF(Py_None);
	py_frontMatBind = NULL;
    }
    if (py_frontMatBind) {
	fmb_array = contiguous_typed_array(py_frontMatBind, PyArray_INT, 1,
					   &expected_bind_dims);
	if (!fmb_array) {
	    Py_DECREF((PyObject *) c_array);
	    Py_DECREF((PyObject *) ind_array);
	    if (tex_array) {
		Py_DECREF((PyObject *) tex_array);
	    }
	    if (n_array) {
		Py_DECREF((PyObject *) n_array);
	    }
	    return NULL;
	}
	fmb = (int *) fmb_array->data;

	for (i = 0; i < 5; i++)
	    frontMatBind[i] = fmb[i];
	hasFrontMatBind = 1;
    }
    /* get the front materials */
    if (py_frontMat == Py_None) {
	Py_DECREF(Py_None);
	py_frontMat = NULL;
    }
    if (py_frontMat) {
	int nd;
	for (i = 0; i < 5; i++) {
	    nd = (i == 4) ? 1 : 2;
	    lstitem = PyList_GetItem(py_frontMat, i);
	    fm_array[i] = contiguous_typed_array(lstitem, PyArray_FLOAT, nd,
						 expected_prop_dims);
	    if (!fm_array[i]) {
		Py_DECREF((PyObject *) c_array);
		Py_DECREF((PyObject *) ind_array);
		if (tex_array) {
		    Py_DECREF((PyObject *) tex_array);
		}
		if (n_array) {
		    Py_DECREF((PyObject *) n_array);
		}
		if (fmb_array) {
		    Py_DECREF((PyObject *) fmb_array);
		}
		for (j = 0; j < i; j++)
		    if (fm_array[j]) {
			Py_DECREF((PyObject *) fm_array[j]);
		    }
		return NULL;
	    }
	    frontMat[i] = (float *) fm_array[i]->data;

	    if (!hasFrontMatBind) {
		if (fm_array[i]->dimensions[0] == c_array->dimensions[0])
		    frontMatBind[i] = PER_VERTEX;
		else if (fm_array[i]->dimensions[0] == ind_array->dimensions[0])
		    frontMatBind[i] = PER_PART;
		else if (fm_array[i]->dimensions[0] == 1)
		    frontMatBind[i] = OVERALL;
		else
		    fm_array[i]->dimensions[0] = NONE;
	    }
	}
    } else {
	for (j = 0; j < 5; j++)
	    fm_array[j] = NULL;
    }


    /* check back material binding parameter */
    if (py_backMatBind == Py_None) {
	Py_DECREF(Py_None);
	py_backMatBind = NULL;
    }
    if (py_backMatBind) {
	bmb_array = contiguous_typed_array(py_backMatBind, PyArray_INT, 1,
					   &expected_bind_dims);
	if (!bmb_array) {
	    Py_DECREF((PyObject *) c_array);
	    Py_DECREF((PyObject *) ind_array);
	    if (tex_array) {
		Py_DECREF((PyObject *) tex_array);
	    }
	    if (n_array) {
		Py_DECREF((PyObject *) n_array);
	    }
	    if (fmb_array) {
		Py_DECREF((PyObject *) fmb_array);
	    }
	    for (j = 0; j < 5; j++)
		if (fm_array[j]) {
		    Py_DECREF((PyObject *) fm_array[j]);
		}
	    return NULL;
	}
	bmb = (int *) bmb_array->data;

	for (i = 0; i < 5; i++)
	    backMatBind[i] = bmb[i];
	hasBackMatBind = 1;
    }
    /* get the back face materials */
    if (py_backMat == Py_None) {
	Py_DECREF(Py_None);
	py_backMat = NULL;
    }
    if (py_backMat) {
	int nd;
	for (i = 0; i < 5; i++) {
	    nd = (i == 4) ? 1 : 2;
	    lstitem = PyList_GetItem(py_backMat, i);
	    bm_array[i] = contiguous_typed_array(lstitem, PyArray_FLOAT, nd,
						 expected_prop_dims);
	    if (!bm_array[i]) {
		Py_DECREF((PyObject *) c_array);
		Py_DECREF((PyObject *) ind_array);
		if (tex_array) {
		    Py_DECREF((PyObject *) tex_array);
		}
		if (n_array) {
		    Py_DECREF((PyObject *) n_array);
		}
		if (fmb_array) {
		    Py_DECREF((PyObject *) fmb_array);
		}
		if (bmb_array) {
		    Py_DECREF((PyObject *) bmb_array);
		}
		for (j = 0; j < 5; j++)
		    if (fm_array[j]) {
			Py_DECREF((PyObject *) fm_array[j]);
		    }
		for (j = 0; j < i; j++)
		    if (bm_array[j]) {
			Py_DECREF((PyObject *) bm_array[j]);
		    }
		return NULL;
	    }
	    backMat[i] = (float *) bm_array[i]->data;

	    if (!hasBackMatBind) {
		if (bm_array[i]->dimensions[0] == c_array->dimensions[0])
		    backMatBind[i] = PER_VERTEX;
		else if (bm_array[i]->dimensions[0] == ind_array->dimensions[0])
		    backMatBind[i] = PER_PART;
		else if (bm_array[i]->dimensions[0] == 1)
		    backMatBind[i] = OVERALL;
		else
		    bm_array[i]->dimensions[0] = NONE;
	    }
	}
    } else {
	for (j = 0; j < 5; j++)
	    bm_array[j] = NULL;
    }

    dpl = glGenLists(1);
    glNewList(dpl, GL_COMPILE_AND_EXECUTE);

    if (!frontAndBack) {
	face = GL_FRONT;
    } else {
	face = GL_FRONT_AND_BACK;
    }

    if (normBinding == NONE) {
	glDisable(GL_LIGHTING);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    } else {
	glEnable(GL_LIGHTING);
	if (normBinding == OVERALL)
	    glNormal3dv(&norm[0]);	/* Overall Normal */
    }

    l = ind_array->dimensions[1];
    if (type == GL_LINES || type == GL_TRIANGLES || type == GL_QUADS) {
	fixed = 1;
	glBegin((GLenum) type);	/* fixed length geom's ==> vertices cannot be picked */
    }
    /* set OVERALL color/material properties
       if (normBinding == NONE) {
       if (py_frontMat) {
       if (frontMatBind[noCol] == OVERALL)
       glColor4fv( &frontMat[noCol][0] );
       }
       } else {
       if (py_frontMat) {
       for (j=0; j<5; j++) {
       if (frontMatBind[j] == OVERALL)
       glMaterialfv( face, propConst[j], &frontMat[j][0] );
       }
       }
       if (py_backMat && !frontAndBack) {
       for (j=0; j<5; j++) {
       if (backMatBind[j] == OVERALL)
       glMaterialfv( GL_BACK, propConst[j], &backMat[j][0] );
       }
       }
       }
     */

    /* initialize color memory */
    isNewColor(NULL);
    isNewMaterial(0, 0, NULL);

    /* loop over faces */
    for (i = 0; i < ind_array->dimensions[0]; i++) {
	if (normBinding == PER_PART)
	    glNormal3dv(&norm[i * 3]);	/* PER_PART */

	/* set PER_PART color/material properties */
	if (normBinding == NONE) {
	    if (py_frontMat) {
		if (frontMatBind[noCol] == PER_PART)
		    if (isNewColor(&frontMat[noCol][i * 4])) {
			glColor4fv(&frontMat[noCol][i * 4]);
		    }
	    }
	} else {
	    if (py_frontMat) {
		for (j = 0; j < 5; j++) {
		    if (frontMatBind[j] == PER_PART)
			if (isNewMaterial(face, j, &frontMat[j][i * 4])) {
			    glMaterialfv((GLenum) face, (GLenum) propConst[j],
					 &frontMat[j][i * 4]);
			}
		}
	    }
	    if (py_backMat && !frontAndBack) {
		for (j = 0; j < 5; j++) {
		    if (backMatBind[j] == PER_PART)
			if (isNewMaterial(GL_BACK, j, &backMat[j][i * 4])) {
			    glMaterialfv(GL_BACK, (GLenum) propConst[j], &backMat[j][i * 4]);
			}
		}
	    }
	}

	if (!fixed) {
	    glPushName(i);
	    glBegin((GLenum) type);
	}
	/* loop over vertices in face */
	for (j = 0; j < l; j++) {
	    v = indices[i * l + j];

	    if (v < 0)
		break;

	    if (v >= c_array->dimensions[0]) {
		sprintf(buf, "Coordinates index %d in face %d out of range", v, j);
		PyErr_SetString(PyExc_ValueError, buf);
		Py_DECREF((PyObject *) c_array);
		Py_DECREF((PyObject *) ind_array);
		if (tex_array) {
		    Py_DECREF((PyObject *) tex_array);
		}
		if (n_array) {
		    Py_DECREF((PyObject *) n_array);
		}
		if (fmb_array) {
		    Py_DECREF((PyObject *) fmb_array);
		}
		if (bmb_array) {
		    Py_DECREF((PyObject *) bmb_array);
		}
		for (i = 0; i < 5; i++) {
		    if (fm_array[i]) {
			Py_DECREF((PyObject *) fm_array[i]);
		    }
		    if (bm_array[i]) {
			Py_DECREF((PyObject *) bm_array[i]);
		    }
		}
		return NULL;
	    }
	    if (normBinding == PER_VERTEX) {
		if (v >= n_array->dimensions[0]) {
		    sprintf(buf, "Normal index %d in face %d out of range", v, j);
		    PyErr_SetString(PyExc_ValueError, buf);
		    Py_DECREF((PyObject *) c_array);
		    Py_DECREF((PyObject *) ind_array);
		    if (tex_array) {
			Py_DECREF((PyObject *) tex_array);
		    }
		    if (n_array) {
			Py_DECREF((PyObject *) n_array);
		    }
		    if (fmb_array) {
			Py_DECREF((PyObject *) fmb_array);
		    }
		    if (bmb_array) {
			Py_DECREF((PyObject *) bmb_array);
		    }
		    for (i = 0; i < 5; i++) {
			if (fm_array[i]) {
			    Py_DECREF((PyObject *) fm_array[i]);
			}
			if (bm_array[i]) {
			    Py_DECREF((PyObject *) bm_array[i]);
			}
		    }
		    return NULL;
		}
		glNormal3dv(&norm[v * 3]);
	    }
	    if (normBinding == NONE) {
		if (py_frontMat) {
		    if (frontMatBind[noCol] == PER_VERTEX)
			if (isNewColor(&frontMat[noCol][v * 4])) {
			    glColor4fv(&frontMat[noCol][v * 4]);
			}
		}
	    } else {
		if (py_frontMat) {
		    for (k = 0; k < 5; k++) {
			if (frontMatBind[k] == PER_VERTEX)
			    if (isNewMaterial(face, k, &frontMat[k][v * 4])) {
				glMaterialfv((GLenum) face, (GLenum) propConst[k],
					     &frontMat[k][v * 4]);
			    }
		    }
		}
		if (py_backMat && !frontAndBack) {
		    for (k = 0; k < 5; k++) {
			if (backMatBind[k] == PER_VERTEX)
			    if (isNewMaterial(GL_BACK, k, &backMat[k][v * 4])) {
				glMaterialfv(GL_BACK, (GLenum) propConst[k],
					     &backMat[k][v * 4]);
			    }
		    }
		}
	    }
	    if (py_texInd) {
		switch (tex_array->dimensions[1]) {
		case 1:
		    glTexCoord1f(texInd[v]);
		    break;
		case 2:
		    glTexCoord2fv(&texInd[v * 2]);
		    break;
		case 3:
		    glTexCoord3fv(&texInd[v * 3]);
		    break;
		case 4:
		    glTexCoord4fv(&texInd[v * 4]);
		    break;
		}
	    }
	    glVertex3dv(&coords[v * 3]);
	}

	if (!fixed) {
	    glEnd();
	    glPopName();
	}
    }

    if (fixed)
	glEnd();

    glEndList();
    Py_DECREF((PyObject *) c_array);
    Py_DECREF((PyObject *) ind_array);
    if (tex_array) {
	Py_DECREF((PyObject *) tex_array);
    }
    if (n_array) {
	Py_DECREF((PyObject *) n_array);
    }
    if (fmb_array) {
	Py_DECREF((PyObject *) fmb_array);
    }
    if (bmb_array) {
	Py_DECREF((PyObject *) bmb_array);
    }
    for (i = 0; i < 5; i++) {
	if (fm_array[i]) {
	    Py_DECREF((PyObject *) fm_array[i]);
	}
	if (bm_array[i]) {
	    Py_DECREF((PyObject *) bm_array[i]);
	}
    }
    return Py_BuildValue("i", dpl);
}

#ifndef _PYMOL_NO_GLUT

#ifndef _PYMOL_OSX
#include <GL/glut.h>
#else
#include <glut.h>
#endif

static PyObject *glSphereSetDSPL(PyObject * self, PyObject * args,
				 PyObject * kw)
{
    static char *argnames[] =
    {"coordinates",
     "frontMaterial", "backMaterial",
     "frontMatBind", "backMatBind",
     "frontAndBack", "noLightCol", "fillMode",
     "slices", "stacks",
     NULL};

    /* int type;*/
      int frontAndBack = 0, noCol = 0, slices = 10, stacks;
    PyObject *py_coords, *py_frontMat = NULL, *py_backMat = NULL, *py_frontMatBind = NULL,
    *py_backMatBind = NULL, *lstitem;

    PyArrayObject *c_array, *fmb_array = NULL, *bmb_array = NULL, *fm_array[5],
    *bm_array[5];

    double *coords;
    float *frontMat[5], *backMat[5];
    int fillMode = -1;
    int expected_coord_dims[2] =
    {-1, 4}, expected_prop_dims[2] =
    {-1, 4};
    int expected_bind_dims = 5, hasFrontMatBind = 0, hasBackMatBind = 0;

    GLuint dpl;
    /* int  k, l, v, w, PER_VERTEX = 11;*/
    int i, j, face, frontMatBind[5], backMatBind[5], *fmb,
    *bmb, NONE = -1, OVERALL = 10, PER_PART = 12, propConst[] =
    {GL_AMBIENT, GL_DIFFUSE, GL_SPECULAR, GL_EMISSION,
     GL_SHININESS};

    if (!PyArg_ParseTupleAndKeywords(args, kw, "O|OOOOiiiii", argnames,
				     &py_coords,
				     &py_frontMat, &py_backMat,
				     &py_frontMatBind, &py_backMatBind,
				     &frontAndBack, &noCol, &fillMode,
				     &slices, &stacks))
	return NULL;

    c_array = contiguous_typed_array(py_coords, PyArray_DOUBLE, 2,
				     expected_coord_dims);
    if (!c_array)
	return NULL;
    coords = (double *) c_array->data;

    /* check front material binding parameter */
    if (py_frontMatBind == Py_None) {
	Py_DECREF(Py_None);
	py_frontMatBind = NULL;
    }
    if (py_frontMatBind) {
	fmb_array = contiguous_typed_array(py_frontMatBind, PyArray_INT, 1,
					   &expected_bind_dims);
	if (!fmb_array) {
	    Py_DECREF((PyObject *) c_array);
	    return NULL;
	}
	fmb = (int *) fmb_array->data;

	for (i = 0; i < 5; i++)
	    frontMatBind[i] = fmb[i];
	hasFrontMatBind = 1;
    }
    /* get the front materials */
    if (py_frontMat == Py_None) {
	Py_DECREF(Py_None);
	py_frontMat = NULL;
    }
    if (py_frontMat) {
	int nd;
	for (i = 0; i < 5; i++) {
	    nd = (i == 4) ? 1 : 2;
	    lstitem = PyList_GetItem(py_frontMat, i);
	    fm_array[i] = contiguous_typed_array(lstitem, PyArray_FLOAT, nd,
						 expected_prop_dims);
	    if (!fm_array[i]) {
		Py_DECREF((PyObject *) c_array);
		if (fmb_array) {
        Py_DECREF((PyObject *) fmb_array);}
		for (j = 0; j < i; j++)
        if (fm_array[j]) {
          Py_DECREF((PyObject *) fm_array[j]);}
		return NULL;
	    }
	    frontMat[i] = (float *) fm_array[i]->data;

	    if (!hasFrontMatBind) {
		if (fm_array[i]->dimensions[0] == c_array->dimensions[0])
		    frontMatBind[i] = PER_PART;
		else if (fm_array[i]->dimensions[0] == 1)
		    frontMatBind[i] = OVERALL;
		else
		    fm_array[i]->dimensions[0] = NONE;
	    }
	}
    }
    /* check back material binding parameter */
    if (py_backMatBind == Py_None) {
	Py_DECREF(Py_None);
	py_backMatBind = NULL;
    }
    if (py_backMatBind) {
	bmb_array = contiguous_typed_array(py_backMatBind, PyArray_INT, 1,
					   &expected_bind_dims);
	if (!bmb_array) {
	    Py_DECREF((PyObject *) c_array);
	    if (fmb_array) {
		Py_DECREF((PyObject *) fmb_array);
       }
	    for (j = 0; j < 5; j++)
         if (fm_array[j]) {
		    Py_DECREF((PyObject *) fm_array[j]);
         }
	    return NULL;
	}
	bmb = (int *) bmb_array->data;

	for (i = 0; i < 5; i++)
	    backMatBind[i] = bmb[i];
	hasBackMatBind = 1;
    }
    /* get the back face materials */
    if (py_backMat == Py_None) {
	Py_DECREF(Py_None);
	py_backMat = NULL;
    }
    if (py_backMat) {
	int nd;
	for (i = 0; i < 5; i++) {
	    nd = (i == 4) ? 1 : 2;
	    lstitem = PyList_GetItem(py_backMat, i);
	    bm_array[i] = contiguous_typed_array(lstitem, PyArray_FLOAT, nd,
						 expected_prop_dims);
	    if (!bm_array[i]) {
		Py_DECREF((PyObject *) c_array);
		if (fmb_array) {
        Py_DECREF((PyObject *) fmb_array);}
		if (bmb_array) {
        Py_DECREF((PyObject *) bmb_array); }
		for (j = 0; j < 5; j++)
        if (fm_array[j]) {
			Py_DECREF((PyObject *) fm_array[j]);
        }
		for (j = 0; j < i; j++)
        if (bm_array[j]) {
			Py_DECREF((PyObject *) bm_array[j]);
        }
		return NULL;
	    }
	    backMat[i] = (float *) bm_array[i]->data;

	    if (!hasBackMatBind) {
		if (bm_array[i]->dimensions[0] == c_array->dimensions[0])
		    backMatBind[i] = PER_PART;
		else if (bm_array[i]->dimensions[0] == 1)
		    backMatBind[i] = OVERALL;
		else
		    bm_array[i]->dimensions[0] = NONE;
	    }
	}
    }
    dpl = glGenLists(1);
    glNewList(dpl, GL_COMPILE_AND_EXECUTE);

    if (!frontAndBack) {
	face = GL_FRONT;
    } else {
	face = GL_FRONT_AND_BACK;
    }

    if (fillMode == GL_LINES) {
	glDisable(GL_LIGHTING);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    } else {
	glEnable(GL_LIGHTING);
    }

    /* set OVERALL color/material properties
       if (fillMode == GL_LINES) {
       if (py_frontMat) {
       if (frontMatBind[noCol] == OVERALL)
       glColor4fv( &frontMat[noCol][0] );
       }
       } else {
       if (py_frontMat) {
       for (j=0; j<5; j++) {
       if (frontMatBind[j] == OVERALL)
       glMaterialfv( face, propConst[j], &frontMat[j][0] );
       }
       }
       if (py_backMat && !frontAndBack) {
       for (j=0; j<5; j++) {
       if (backMatBind[j] == OVERALL)
       glMaterialfv( GL_BACK, propConst[j], &backMat[j][0] );
       }
       }
       }
     */

    /* initialize color memory */
    isNewColor(NULL);
    isNewMaterial(0, 0, NULL);

    /* loop over spheres */
    for (i = 0; i < c_array->dimensions[0]; i++) {

	/* set PER_PART (e.g. per sphere) color/material properties */
	if (fillMode == GL_LINES) {
	    if (py_frontMat) {
		if (frontMatBind[noCol] == PER_PART)
		    if (isNewColor(&frontMat[noCol][i * 4])) {
			glColor4fv(&frontMat[noCol][i * 4]);
		    }
	    }
	} else {
	    if (py_frontMat) {
		for (j = 0; j < 5; j++) {
		    if (frontMatBind[j] == PER_PART) {
			if (isNewMaterial(face, j, &frontMat[j][i * 4])) {
			    glMaterialfv((GLenum) face, (GLenum) propConst[j],
					 &frontMat[j][i * 4]);
			}
		    }
		}
	    }
	    if (py_backMat && !frontAndBack) {
		for (j = 0; j < 5; j++) {
		    if (backMatBind[noCol] == PER_PART) {
			if (isNewMaterial(GL_BACK, j, &backMat[j][i * 4])) {
			    glMaterialfv(GL_BACK, (GLenum) propConst[j], &backMat[j][i * 4]);
			}
		    }
		}
	    }
	}

	glPushName(i);
	glPushMatrix();
	glTranslated(coords[i * 4], coords[i * 4 + 1], coords[i * 4 + 2]);
	glutSolidSphere(coords[i * 4 + 3], slices, stacks);
	glPopMatrix();
	glPopName();
    }
    glEndList();

    Py_DECREF((PyObject *) c_array);
    if (fmb_array) {
      Py_DECREF((PyObject *) fmb_array); }
    if (bmb_array) {
      Py_DECREF((PyObject *) bmb_array); }
    for (j = 0; j < 5; j++) {
      if (fm_array[j]) {
        Py_DECREF((PyObject *) fm_array[j]); }
      if (bm_array[j]) { 
        Py_DECREF((PyObject *) bm_array[j]); }
    }
    return Py_BuildValue("i", dpl);
}
#endif

/****************************************************************
  TRACKBALL Object
****************************************************************/

/*
 * Local function for the trackball
 */

static void track_vcopy(const float *v1, float *v2)
{
    register int i;
    for (i = 0; i < 3; i++)
	v2[i] = v1[i];
}

static void track_vcross(const float *v1, const float *v2, float *cross)
{
    float temp[3];

    temp[0] = (v1[1] * v2[2]) - (v1[2] * v2[1]);
    temp[1] = (v1[2] * v2[0]) - (v1[0] * v2[2]);
    temp[2] = (v1[0] * v2[1]) - (v1[1] * v2[0]);
    track_vcopy(temp, cross);
}

static float track_vlength(const float *v)
{
    return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

static void track_vscale(float *v, float div)
{
    v[0] *= div;
    v[1] *= div;
    v[2] *= div;
}

static void track_vnormal(float *v)
{
    track_vscale(v, 1.0 / track_vlength(v));
}

/*
 *  Given an axis and angle, compute quaternion.
 */
static void track_axis_to_quat(float a[3], float phi, float q[4])
{
    track_vnormal(a);
    track_vcopy(a, q);
    track_vscale(q, sin(phi / 2.0));
    q[3] = cos(phi / 2.0);
}

/*
 * Project an x,y pair onto a sphere of radius r OR a hyperbolic sheet
 * if we are away from the center of the sphere.
 */
static float track_project_to_sphere(float r, float x, float y)
{
    float d, t, z;

    d = sqrt(x * x + y * y);
    if (d < r * 0.70710678118654752440) {	/* Inside sphere */
	z = sqrt(r * r - d * d);
    } else {			/* On hyperbola */
	t = r / 1.41421356237309504880;
	z = t * t / d;
    }
    return z;
}

/*
 * Simulate a track-ball.  Project the points onto the virtual
 * trackball, then figure out the axis of rotation, which is the cross
 * product of P1 P2 and O P1 (O is the center of the ball, 0,0,0)
 * Note:  This is a deformed trackball-- is a trackball in the center,
 * but is deformed into a hyperbolic sheet of rotation away from the
 * center.  This particular function was chosen after trying out
 * several variations.
 *
 * p1x, p1y: last cursor position  (assumed in range -1.0 ... 1.0)
 * p2x, p2y: current cursor position  (assumed in range -1.0 ... 1.0)
 * 
 * result: q, quaternion describing the rotation
 */
static void trackball(float q[4], float p1x, float p1y, float p2x, float p2y,
		      float size)
{
    int i;
    float a[3];			/* Axis of rotation */
    float phi;			/* how much to rotate about axis */
    float p1[3], p2[3], d[3];
    float t;

    if (p1x == p2x && p1y == p2y) {	/* Zero rotation */
	q[0] = q[1] = q[2] = 0.0;
	q[3] = 1.0;
	return;
    }
    /*
     * First, figure out z-coordinates for projection of P1 and P2 to
     * deformed sphere
     */
    p1[0] = p1x;
    p1[1] = p1y;
    p1[2] = track_project_to_sphere(size, p1x, p1y);
    p2[0] = p2x;
    p2[1] = p2y;
    p2[2] = track_project_to_sphere(size, p2x, p2y);

    /*
     *  Now, we want the cross product of P1 and P2
     */
    track_vcross(p2, p1, a);

    /*
     *  Figure out how much to rotate around that axis.
     */
    for (i = 0; i < 3; i++)
	d[i] = p1[i] - p2[i];

    t = track_vlength(d) / (2.0 * size);

    /*
     * Avoid problems with out-of-control values...
     */
    if (t > 1.0)
	t = 1.0;
    if (t < -1.0)
	t = -1.0;
    phi = 2.0 * asin(t);

    track_axis_to_quat(a, phi, q);
}

/*
 * Build a rotation matrix, given a quaternion rotation.
 *
 */
static void track_build_rotmatrix(float m[4][4], float q[4])
{
    m[0][0] = 1.0 - 2.0 * (q[1] * q[1] + q[2] * q[2]);
    m[0][1] = 2.0 * (q[0] * q[1] - q[2] * q[3]);
    m[0][2] = 2.0 * (q[2] * q[0] + q[1] * q[3]);
    m[0][3] = 0.0;

    m[1][0] = 2.0 * (q[0] * q[1] + q[2] * q[3]);
    m[1][1] = 1.0 - 2.0 * (q[2] * q[2] + q[0] * q[0]);
    m[1][2] = 2.0 * (q[1] * q[2] - q[0] * q[3]);
    m[1][3] = 0.0;

    m[2][0] = 2.0 * (q[2] * q[0] - q[1] * q[3]);
    m[2][1] = 2.0 * (q[1] * q[2] + q[0] * q[3]);
    m[2][2] = 1.0 - 2.0 * (q[1] * q[1] + q[0] * q[0]);
    m[2][3] = 0.0;

    m[3][0] = 0.0;
    m[3][1] = 0.0;
    m[3][2] = 0.0;
    m[3][3] = 1.0;
}

/*****************************************************************/

#define Py_Try(BOOLEAN) {if(!(BOOLEAN)) return NULL;}

#define PyObjtrackball_Check(op) ((op)->ob_type == &PyObjtrackball_type)

typedef struct {
    PyObject_HEAD
    float trackballsize;
    float scale;
    float quat[4];
    float matrix[4][4];
    int renormcount;
} PyObjtrackball;

staticforward PyTypeObject PyObjtrackball_type;

static PyObject *Pytrackball(PyObject * self, PyObject * args)
{
    float p1x, p1y, p2x, p2y;
    int width, height, mat = 0;
    PyObjtrackball *t = (PyObjtrackball *) self;

    if (!PyArg_ParseTuple(args, "ffffii|i", &p1x, &p1y, &p2x, &p2y,
			  &width, &height, &mat))
	return NULL;

    trackball(t->quat,
	      (t->scale * p1x - width) / width,		/* assumed to be -1.0 .... 1.0 */
	      (height - t->scale * p1y) / height,
	      (t->scale * p2x - width) / width,
	      (height - t->scale * p2y) / height,
	      t->trackballsize);

    if (mat)
	track_build_rotmatrix(t->matrix, t->quat);

    Py_INCREF(Py_None);
    return Py_None;
}

static struct PyMethodDef PyObjtrackball_methods[] =
{
    {"update", Pytrackball, 1},
    {NULL, NULL}
};

static PyObject *PyObjtrackball_repr(PyObjtrackball * self)
{
    return Py_BuildValue("s", "Trackball Object");
}

static PyObjtrackball *Newtrackball(float size, float scale, int renorm)
{
    PyObjtrackball *self;
    int i, j;

    Py_Try((self = PyObject_NEW(PyObjtrackball, &PyObjtrackball_type)));

/******Initialize your structure values here******/
    self->trackballsize = size;
    self->scale = scale;
    self->renormcount = renorm;
    for (i = 0; i < 4; i++) {
	self->quat[i] = 0.0;
	for (j = 0; j < 4; j++) {
	    self->matrix[i][j] = 0.0;
	}
	self->matrix[i][i] = 1.0;
    }
/*************************************************/

    return self;
}
static void PyObjtrackball_dealloc(PyObjtrackball * self)
{
    PyMem_DEL(self);
}

static int PyObjtrackball_print(PyObjtrackball * self)
{
    printf("  size  : %f\n", self->trackballsize);
    printf("  scale : %f\n", self->scale);
    printf("  renorm: %i\n", self->renormcount);
    printf("  quat  : %6.3f %6.3f %6.3f %6.3f\n",
	   self->quat[0], self->quat[1], self->quat[2], self->quat[3]);
    printf("  mat   : %6.3f %6.3f %6.3f %6.3f\n",
	   self->matrix[0][0], self->matrix[0][1],
	   self->matrix[0][2], self->matrix[0][3]);
    printf("          %6.3f %6.3f %6.3f %6.3f\n",
	   self->matrix[1][0], self->matrix[1][1],
	   self->matrix[1][2], self->matrix[1][3]);
    printf("          %6.3f %6.3f %6.3f %6.3f\n",
	   self->matrix[2][0], self->matrix[2][1],
	   self->matrix[2][2], self->matrix[2][3]);
    printf("          %6.3f %6.3f %6.3f %6.3f\n",
	   self->matrix[3][0], self->matrix[3][1],
	   self->matrix[3][2], self->matrix[3][3]);

    return 0;
}

/*
   return an Numeric 1D array of 'len' floats
 */
static PyArrayObject *track_array_vector_float(float *data, int len)
{
    PyArrayObject *vector;

    vector = (PyArrayObject *) PyArray_FromDims(1, &len, PyArray_FLOAT);
    if (!vector) {
	PyErr_SetString(PyExc_RuntimeError,
			"Failed to allocate memory for vector");
	return NULL;
    }
    memcpy(vector->data, data, len * sizeof(float));
    return vector;
}

static PyObject *PyObjtrackball_getattr(PyObjtrackball * self, char *name)
{
    if (strcmp(name, "size") == 0)
	return Py_BuildValue("f", self->trackballsize);
    if (strcmp(name, "scale") == 0)
	return Py_BuildValue("f", self->scale);
    else if (strcmp(name, "quat") == 0)
	return (PyObject *) (track_array_vector_float(self->quat, 4));
    else if (strcmp(name, "mat") == 0)
	return (PyObject *) (track_array_vector_float((float *) (self->matrix), 16));
    else if (strcmp(name, "renorm") == 0)
	return Py_BuildValue("i", self->renormcount);

    return Py_FindMethod(PyObjtrackball_methods,
			 (PyObject *) self,
			 name);
}

static int PyObjtrackball_setattr(PyObjtrackball * self, char *name, PyObject * v)
{
    if (strcmp(name, "size") == 0) {
      PyArg_Parse(v, "f", &self->trackballsize);
	return 0;
    } else if (strcmp(name, "scale") == 0) {
      PyArg_Parse(v, "f", &self->scale);
	return 0;
    } else if (strcmp(name, "renom") == 0) {
      PyArg_Parse(v, "i", &self->renormcount);
	return 0;
    }
    PyErr_SetString(PyExc_ValueError, "Sorry, bad or ReadOnly data member");
    return 1;
}

static PyObject *Create_trackball(PyObject * self, PyObject * args, PyObject * kw)
{
    static char *argnames[] =
    {"size", "scale", "renorm", NULL};
    PyObjtrackball *result;
    int renorm = 97;
    float size = 0.8, scale = 2.0;

    if (!PyArg_ParseTupleAndKeywords(args, kw, "|ffi", argnames,
				     &size, &scale, &renorm))
	return NULL;

    result = Newtrackball(size, scale, renorm);

    if (!result) {
	PyErr_SetString(PyExc_RuntimeError, "Failed to allocate memory");
	return NULL;
    }
    return (PyObject *) result;
}


static PyTypeObject PyObjtrackball_type =
{
#ifdef MS_WIN32
    PyObject_HEAD_INIT(NULL)
#else
    PyObject_HEAD_INIT(&PyType_Type)
#endif
    0,				/* Object size              */
    "trackball",
    sizeof(PyObjtrackball),
    0,				/* Item size                */
    (destructor) PyObjtrackball_dealloc,
    (printfunc) PyObjtrackball_print,
    (getattrfunc) PyObjtrackball_getattr,
    (setattrfunc) PyObjtrackball_setattr,
    (cmpfunc) 0,		/* Comparing method         */
    (reprfunc) PyObjtrackball_repr,
    0,				/* As number                */
    0,				/* As sequence              */
    0,				/* As mapping               */
    (hashfunc) 0,		/* Hash function            */
    (ternaryfunc) 0,		/* Ternary function (call)  */
    (reprfunc) 0,		/* Unknown                  */
    0L, 0L, 0L, 0L,		/* Free space               */
    0L				/* Documentation            */
};

/* END MS */

static PyObject *gl_TrianglesWithNormals(PyObject * self, PyObject * args)
{
    int i, size;
    double d, d1[3], d2[3], norm[3], *v;
    PyObject *vop;
    PyArrayObject *mp;
    char *data;
    TRY(PyArg_ParseTuple(args, "O", &vop));
    TRY(mp =
	(PyArrayObject *) PyArray_ContiguousFromObject(vop, PyArray_DOUBLE, 1, 0));
    size = PyArray_Size((PyObject *) mp);
    if ((size % 9) != 0) {
	PyErr_SetString(PyExc_ValueError, "matrix length sould be divisible by 9");
	return NULL;
    }
    glBegin(GL_TRIANGLES);
    for (data = mp->data, i = 0; i < size; i += 9, data += 9 * mp->descr->elsize) {
	v = (double *) data;
	d1[0] = v[0] - v[3];
	d2[0] = v[3] - v[6];
	d1[1] = v[1] - v[4];
	d2[1] = v[4] - v[7];
	d1[2] = v[2] - v[5];
	d2[2] = v[5] - v[8];
	norm[0] = d1[1] * d2[2] - d1[2] * d2[1];
	norm[1] = d1[2] * d2[0] - d1[0] * d2[2];
	norm[2] = d1[0] * d2[1] - d1[1] * d2[0];
	d = sqrt(norm[0] * norm[0] + norm[1] * norm[1] + norm[2] * norm[2]);
	norm[0] /= d;
	norm[1] /= d;
	norm[2] /= d;
	glNormal3dv(norm);
	glVertex3dv(v);
	glVertex3dv(v + 3);
	glVertex3dv(v + 6);
    }
    glEnd();
    Py_DECREF(mp);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *gl_Triangles(PyObject * self, PyObject * args)
{
    int i, size;
    PyObject *vop;
    PyArrayObject *mp;
    char *data;
    double *v;
    TRY(PyArg_ParseTuple(args, "O", &vop));
    TRY(mp =
	(PyArrayObject *) PyArray_ContiguousFromObject(vop, PyArray_DOUBLE, 1, 0));
    size = PyArray_Size((PyObject *) mp);
    if ((size % 9) != 0) {
	PyErr_SetString(PyExc_ValueError, "matrix length sould be divisible by 9");
	return NULL;
    }
    glBegin(GL_TRIANGLES);
    for (data = mp->data, i = 0; i < size; i += 9, data += 9 * mp->descr->elsize) {
	v = (double *) data;
	glVertex3dv(v);
	glVertex3dv(v + 3);
	glVertex3dv(v + 6);
    }
    glEnd();
    Py_DECREF(mp);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *gl_Lines(PyObject * self, PyObject * args)
{
    int i, size;
    PyObject *vop;
    PyArrayObject *mp;
    char *data;
    TRY(PyArg_ParseTuple(args, "O", &vop));
    TRY(mp = (PyArrayObject *) PyArray_ContiguousFromObject(vop, PyArray_DOUBLE, 1, 0));
    size = PyArray_Size((PyObject *) mp);
    if ((size % 3) != 0) {
	PyErr_SetString(PyExc_ValueError, "matrix length sould be divisible by 3");
	return NULL;
    }
    glBegin(GL_LINES);
    for (data = mp->data, i = 0; i < size; i += 3, data += 3 * mp->descr->elsize)
	glVertex3dv((double *) data);
    glEnd();
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *gl_Points(PyObject * self, PyObject * args)
{
    int i, size;
    PyObject *vop;
    PyArrayObject *mp;
    char *data;
    TRY(PyArg_ParseTuple(args, "O", &vop));
    TRY(mp = (PyArrayObject *) PyArray_ContiguousFromObject(vop, PyArray_DOUBLE, 1, 0));
    size = PyArray_Size((PyObject *) mp);
    if ((size % 3) != 0) {
	PyErr_SetString(PyExc_ValueError, "matrix length sould be divisible by 3");
	return NULL;
    }
    glBegin(GL_POINTS);
    for (data = mp->data, i = 0; i < size; i += 3, data += 3 * mp->descr->elsize)
	glVertex3dv((double *) data);
    glEnd();
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *gl_Tetrahedra(PyObject * self, PyObject * args)
{
    int i, size;
    PyObject *vop;
    PyArrayObject *mp;
    char *data;
    double *v;
    TRY(PyArg_ParseTuple(args, "O", &vop));
    TRY(mp = (PyArrayObject *) PyArray_ContiguousFromObject(vop, PyArray_DOUBLE, 1, 0));
    size = PyArray_Size((PyObject *) mp);
    if ((size % 12) != 0) {
	PyErr_SetString(PyExc_ValueError, "matrix length sould be divisible by 12");
	return NULL;
    }
    glBegin(GL_LINE_LOOP);
    for (data = mp->data, i = 0; i < size; i += 12, data += 12 * mp->descr->elsize) {
	v = (double *) data;
	glVertex3dv(v);
	glVertex3dv(v + 3);
	glVertex3dv(v + 6);
	glVertex3dv(v);
	glVertex3dv(v + 9);
	glVertex3dv(v + 3);
	glVertex3dv(v + 6);
	glVertex3dv(v + 9);
    }
    glEnd();
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *gl_Vertex(PyObject * self, PyObject * args)
{
    PyObject *vop;
    PyArrayObject *ap;
    double *dp, *max_dp;

    TRY(PyArg_ParseTuple(args, "O", &vop));

    TRY(ap = (PyArrayObject *) PyArray_ContiguousFromObject(vop, PyArray_DOUBLE, 1, 2));

    max_dp = dp = (double *) ap->data;
    max_dp += PyArray_Size((PyObject *) ap);

    switch (ap->dimensions[1]) {
    case 2:
	while (dp < max_dp) {
	    glVertex2dv(dp);
	    dp += 2;
	}
	break;
    case 3:
	while (dp < max_dp) {
	    glVertex3dv(dp);
	    dp += 3;
	}
	break;
    case 4:
	while (dp < max_dp) {
	    glVertex4dv(dp);
	    dp += 4;
	}
	break;
    default:
	ASSERT(0, "1-4d vertices required");
    }

    Py_DECREF(ap);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *gl_CallLists(PyObject * self, PyObject * args)
{
    PyObject *op;
    PyArrayObject *ap;
    int type;

    TRY(PyArg_ParseTuple(args, "O", &op));
    if (PyArray_Check(op)) {
	ap = (PyArrayObject *) op;
	if (ap->nd != 1) {
	    PyErr_SetString(gl_Error, "calllists: array must be 1d");
	    return NULL;
	}
	Py_INCREF(ap);
    } else {
	TRY(ap = (PyArrayObject *) PyArray_ContiguousFromObject(op, PyArray_INT, 1, 1));
    }

    type = typecode2gltype[ap->descr->type_num];
    if (type == -1) {
	PyErr_SetString(gl_Error, "can't callists on this type of array");
	Py_DECREF(ap);
	return NULL;
    }
    glCallLists(ap->dimensions[0], (GLenum) type, ap->data);

    Py_DECREF(ap);
    Py_INCREF(Py_None);
    return Py_None;
}

#endif


/* ###################### TOGL STUFF ###################################### */

static PyObject *gl_TranslateScene(PyObject * self, PyObject * args)
{
    GLint x, y, mousex, mousey;
    double s;
    GLdouble mat[16];
    if (!PyArg_ParseTuple(args, "diiii", &s, &x, &y, &mousex, &mousey))
	return NULL;
    glMatrixMode(GL_MODELVIEW);
    glGetDoublev(GL_MODELVIEW_MATRIX, mat);
    glLoadIdentity();
    glTranslatef((float)(s * (x - mousex)),(float)( s * (mousey - y)), 0.0F);
    glMultMatrixd(mat);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *gl_RotateScene(PyObject * self, PyObject * args)
{
    GLint x, y, mousex, mousey;
    GLdouble s, xcenter, ycenter, zcenter, mat[16];
    if (!PyArg_ParseTuple(args, "ddddiiii", &s, &xcenter, &ycenter, &zcenter,
			  &x, &y, &mousex, &mousey))
	return NULL;
    glMatrixMode(GL_MODELVIEW);
    glGetDoublev(GL_MODELVIEW_MATRIX, mat);
    glLoadIdentity();
    glTranslatef((float)xcenter,(float) ycenter,(float) zcenter);
    glRotatef((float)(s * (y - mousey)), 1.0F, 0.0F, 0.0F);
    glRotatef((float)(s * (x - mousex)), 0.0F, 1.0F, 0.0F);
    glTranslatef((float)-xcenter,(float) -ycenter,(float) -zcenter);
    glMultMatrixd(mat);
    Py_INCREF(Py_None);
    return Py_None;
}


#define v3op(a,op,b)      (a[0] op b[0], a[1] op b[1], a[2] op b[2])
#define v3eqn(a,aop,b,o,c)(a[0] aop b[0] o c[0], a[1] aop b[1] o c[1], \
                           a[2] aop b[2] o c[2])
#define v3dot(a,b)        ((a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2]))
#define v3lengthsq(v)     v3dot(v, v)
#define VSQ(a)            ((a) * (a))
#define v3distsq(a,b)     (VSQ(a[0] - b[0]) + VSQ(a[1] - b[1]) + VSQ(a[2] - b[2]))

/*
 * Square of perpendicular distance of point from line.
 * gl.DistFromLine(x, p1, p2)
 * p1 and p2 are the end points of the line (for instance
 * returned from gl.UnProject.
 * Should transform this to Python documentation!
 */

static PyObject *gl_DistFromLine(PyObject * self, PyObject * args)
{
    PyObject *ox, *op1, *op2, *ret;
    GLdouble *x, *p1, *p2;
    int nx, np1, np2;
    GLdouble p1top2[3], p2top1[3], p1tox[3], p2tox[3];
    GLdouble dsq, dot;

    if (!PyArg_ParseTuple(args, "OOO", &ox, &op1, &op2))
	return NULL;

    TRY(PyArray_AsDoubleArray(&ox, &x, &nx));
    TRY(PyArray_AsDoubleArray(&op1, &p1, &np1));
    TRY(PyArray_AsDoubleArray(&op2, &p2, &np2));

    if (nx != 3 || np1 != 3 || np2 != 3) {
	PyErr_SetString(gl_Error, "All arguments should contain 3 items!");
	PyArray_ClearMemory(ox, x);
	PyArray_ClearMemory(op1, p1);
	PyArray_ClearMemory(op2, p2);
	return NULL;
    }
/* check if point lies off either end of line segment */

    v3eqn(p1top2, =, p2, -, p1);
    v3eqn(p1tox, =, x, -, p1);

    if (v3dot(p1tox, p1top2) < 0.0) {
	PyArray_ClearMemory(ox, x);
	PyArray_ClearMemory(op1, p1);
	PyArray_ClearMemory(op2, p2);
	return PyFloat_FromDouble(HUGE_VAL);	/* was HUGE daa */
    }
    v3op(p2top1, =, -p1top2);
    v3eqn(p2tox, =, x, -, p2);

    if (v3dot(p2tox, p2top1) < 0.0) {
	PyArray_ClearMemory(ox, x);
	PyArray_ClearMemory(op1, p1);
	PyArray_ClearMemory(op2, p2);
	return PyFloat_FromDouble(HUGE_VAL);	/* was HUGE daa */
    }
/* from dot product */

    dot = v3dot(p1top2, p1tox);

    dsq = dot * dot / v3lengthsq(p1top2);

/* from pythagoras */

    ret = PyFloat_FromDouble(v3distsq(p1, x) - dsq);
    PyArray_ClearMemory(ox, x);
    PyArray_ClearMemory(op1, p1);
    PyArray_ClearMemory(op2, p2);
    return ret;
}



/* ########################## END OF TOGL STUFF ########################### */

static PyMethodDef glutil_methods[] =
{

#ifdef LIBTIFF
    {"glSaveTiff", gl_SaveTiff, 1},
#endif
    {"glSavePPM", gl_SavePPM, 1},
    {"glSaveEPS", gl_SaveEPS, 1},

/* MS Dec 12 '97 */
#ifdef NUMERIC
    {"glCleanRotMat", gl_CleanRotMat, 1},
    {"glTriangleNormals", gl_TriangleNormals, 1},
    {"glIndexedGeomDSPL", (PyCFunction) gl_indexedGeomDSPL, 3},
#ifndef _PYMOL_NO_GLUT
    {"glSphereSetDSPL", (PyCFunction) glSphereSetDSPL, 3},
#endif
    {"glTrackball", (PyCFunction) Create_trackball, 3},
/* END MS */

    {"glTrianglesWithNormals", gl_TrianglesWithNormals, 1},
    {"glTriangles", gl_Triangles, 1},
    {"glLines", gl_Lines, 1},
    {"glPoints", gl_Points, 1},
    {"glTetrahedra", gl_Tetrahedra, 1},
    {"glVertex", gl_Vertex, 1},
    {"glCallLists", gl_CallLists, 1},
#endif				/* Not NUMERIC */

    {"glColorVertex2", gl_ColorVertex2d, 1},
    {"glColorVertex2d", gl_ColorVertex2d, 1},

    {"glTranslateScene", gl_TranslateScene, 1},
    {"glRotateScene", gl_RotateScene, 1},
    {"glDistFromLine", gl_DistFromLine, 1},
    {"glSelectWithCallback", gl_SelectWithCallback, 1, gl_SelectWithCallback__doc__},

    {NULL, NULL}
};

static char openglutil_module_documentation[] =
"OpenGl relevant functions:\n\
   glSaveTiff (if LIBTIFF)\n\
   glSavePPM\n\
   glSaveEPS\n\
   glGetError\n\
 #ifdef NUMERIC\n\
   glSelectBuffer\n\
   glCleanRotMat\n\
   glTriangleNormals\n\
   glindexedGeomDSPL\n\
   gltrackball\n\
   glTrianglesWithNormals\n\
   glTriangles\n\
   glLines\n\
   glPoints\n\
   glTetrahedra\n\
   glVertex\n\
   glCallLists\n\
 #endif\n\
   glColorVertex2\n\
   glColorVertex2d\n\
   glTranslateScene\n\
   glRotateScene\n\
   glDistFromLine\n\
   glSelectWithCallback";

#ifdef NUMERIC
DL_EXPORT(void) initopenglutil_num(void);
DL_EXPORT(void) initopenglutil_num(void)
#else
DL_EXPORT(void) initopenglutil(void);
DL_EXPORT(void) initopenglutil(void)
#endif
{
    PyObject *m, *d;
    PyObject *nl, *gl;

#ifdef NUMERIC
#ifdef MS_WIN32
    PyObjtrackball_type.ob_type = &PyType_Type;
#endif
    m = Py_InitModule4("openglutil_num",
#else
    m = Py_InitModule4("openglutil",
#endif
		       glutil_methods,
		       openglutil_module_documentation,
		       (PyObject *) NULL,
		       PYTHON_API_VERSION);

    d = PyModule_GetDict(m);
#ifdef NUMERIC
#ifdef import_array		/* Konrad Hinsen's version */
    import_array();
#endif
    gl_Error = Py_BuildValue("s", "openglutil_num.error");
#else
    gl_Error = Py_BuildValue("s", "openglutil.error");
#endif
    PyDict_SetItemString(d, "error", gl_Error);
#ifdef NUMERIC
    nl = PyInt_FromLong(1L);
#else
    nl = PyInt_FromLong(0L);
#endif
    PyDict_SetItemString(d, "_numeric", nl);
    Py_DECREF(nl);
#ifndef PYMOL_NO_GLUT
    gl = PyInt_FromLong(1L);
#else
    gl = PyInt_FromLong(0L);
#endif
    PyDict_SetItemString(d, "_glut", gl);
    Py_DECREF(gl);
    if (PyErr_Occurred())
	Py_FatalError("can't initialize module openglutil");
}
