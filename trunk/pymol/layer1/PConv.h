
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
#ifndef _H_PConv
#define _H_PConv

#include"os_python.h"

#ifndef _PYMOL_NOPY
#include"Base.h"
#include"OVLexicon.h"


/* Convenient conversion routines for C<->Python data interchange
   
   Note that all of these routines assume that we have the global
   interpreter lock - blocking all other threads.
   
   There are three ways to get it:
   
   - call PBlock() [followe by PUnblock() when done]
   
   - call PBlockAndUnlockAPI - [followed by PLockAPIAndUnblock() when
   done]
   
   - or in response to a call to the PM API, you will have the main
   Python thread by default.  [Note that within an
   APIEntry(),APIExit() block the lock is released, so these
   functions should be called outside of that block].

*/


/* == error-checking routines: true = success, false = failure. */


/* NOTE: the string routines will write strings up to the specified
 * length, PLUS a NULL...so watch out for array overruns */

int PConvAttrToStrMaxLen(PyObject * obj, char *attr, char *str, ov_size ll);

int PConvPyListToExtent(PyObject * obj, float *mn, float *mx);

int PConvAttrToFloatArrayInPlace(PyObject * obj, char *attr, float *ff, ov_size ll);
int PConvAttrToIntArrayInPlace(PyObject * obj, char *attr, int *ff, ov_size ll);
int PConvAttrToPtr(PyObject * obj, char *name, void **cobj);

int PConvCObjectToPtr(PyObject * obj, void **ptr);
int PConvPyListToStrVLAList(PyObject * obj, char **vla, int *n_str);

int PConvPyListToStringVLA(PyObject * obj, char **vla_ptr);
int PConvPyListToIntVLA(PyObject * obj, int **f);
int PConvPyStrToStr(PyObject * obj, char *ptr, int l);
int PConvPyStrToStrPtr(PyObject * obj, char **ptr);
int PConvPyStrToLexRef(PyObject * obj, OVLexicon * lex, int *lex_ref);
int PConvPyFloatToFloat(PyObject * obj, float *ptr);
int PConvPyIntToChar(PyObject * obj, char *ptr);
int PConvPyIntToInt(PyObject * obj, int *ptr);
int PConvPyListToLabPosVLA(PyObject * obj, LabPosType ** vla_ptr);


/* Jenarix conventions -- returns before args */

ov_status PConvPyTupleToIntVLA(int **result, PyObject * tuple);
ov_status PConvPyTupleToFloatVLA(float **result, PyObject * tuple);


/* === end === */


/* categories below... */

PyObject *PConvFloatVLAToPyList(float *vla);
PyObject *PConvFloatVLAToPyTuple(float *vla);
PyObject *PConvIntVLAToPyList(int *vla);
PyObject *PConvIntVLAToPyTuple(int *vla);
PyObject *PConvIntArrayToPyList(int *f, int l);
PyObject *PConvSIntArrayToPyList(short int *f, int l);
PyObject *PConvSCharArrayToPyList(signed char *f, int l);
PyObject *PConvLabPosVLAToPyList(LabPosType * vla, int l);


/* WARNING: the returned PyObject is unowned - it is intended for use
 * only for efficient detection of changes to dictionary values
 * following evaluation of some expression in the context of the
 * dictionary PAlter, PAlterState, etc. */
PyObject *PConvFloatToPyDictItem(PyObject * dict, char *key, float f);
PyObject *PConvStringToPyDictItem(PyObject * dict, char *key, char *f);
PyObject *PConvIntToPyDictItem(PyObject * dict, char *key, int i);

/* end WARNING */

void PConvFloat3ToPyObjAttr(PyObject * obj, char *attr, float *v);
void PConvFloatToPyObjAttr(PyObject * obj, char *attr, float f);
void PConvIntToPyObjAttr(PyObject * obj, char *attr, int i);
void PConvInt2ToPyObjAttr(PyObject * obj, char *attr, int *v);
void PConvStringToPyObjAttr(PyObject * obj, char *attr, char *f);

int PConvPyObjectToFloat(PyObject * object, float *value);
int PConvPyObjectToInt(PyObject * object, int *value);
int PConvPyObjectToChar(PyObject * object, char *value);


/* NOTE: the string routines will write strings up to the specified
 * length, PLUS a NULL...so watch out for array overruns */

int PConvPyObjectToStrMaxLen(PyObject * object, char *value, int ln);
int PConvPyObjectToStrMaxClean(PyObject * object, char *value, int ln);

PyObject *PConvStringListToPyList(int l, char **str);
PyObject *PConvStringVLAToPyList(char *str);

void PConv44PyListTo44f(PyObject * src, float *dest);   /* note loss of precision */

int PConvPyListToFloatVLA(PyObject * obj, float **f);
int PConvPyListToFloatVLANoneOkay(PyObject * obj, float **f);
int PConvPyList3ToFloatVLA(PyObject * obj, float **f);
int PConvPyListToFloatArray(PyObject * obj, float **f);
int PConvPyListToDoubleArray(PyObject * obj, double **f);
int PConvPyListToFloatArrayInPlace(PyObject * obj, float *ff, ov_size ll);
int PConvPyListToFloatArrayInPlaceAutoZero(PyObject * obj, float *ii, ov_size ll);

int PConvPyListToDoubleArrayInPlace(PyObject * obj, double *ff, ov_size ll);

PyObject *PConvFloatArrayToPyList(float *f, int l);
PyObject *PConvFloatArrayToPyListNullOkay(float *f, int l);
PyObject *PConvDoubleArrayToPyList(double *f, int l);

int PConvPyListToIntArray(PyObject * obj, int **f);
int PConvPyListToIntArrayInPlace(PyObject * obj, int *ff, ov_size ll);
int PConvPyListToIntArrayInPlaceAutoZero(PyObject * obj, int *ii, ov_size ll);

int PConvPyListToSIntArrayInPlaceAutoZero(PyObject * obj, short int *ii, ov_size ll);
int PConvPyListToSCharArrayInPlaceAutoZero(PyObject * obj, signed char *ii, ov_size ll);

PyObject *PConv3DIntArrayTo3DPyList(int ***array, int *dim);

PyObject *PConvAutoNone(PyObject * result);     /* automatically own Py_None */

#endif
#endif
