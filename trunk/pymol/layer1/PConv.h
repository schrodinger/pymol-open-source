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

#include<Python.h>


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

int PConvAttrToStrMaxLen(PyObject *obj,char *attr,char *str,int ll);

int PConvPyListToExtent(PyObject *obj,float *mn,float *mx);

int PConvAttrToFloatArrayInPlace(PyObject *obj,char *attr,float *ff,int ll);
int PConvAttrToIntArrayInPlace(PyObject *obj,char *attr,int *ff,int ll);
int PConvAttrToPtr(PyObject *obj,char *name,void **cobj);

int PConvCObjectToPtr(PyObject *obj,void **ptr);

int PConvPyListToStringVLA(PyObject *obj,char **vla_ptr);

/* === end === */

/* === optimized, non-error checking routines === */

/* === end === */

/* categories below... */

PyObject *PConvFloatVLAToPyList(float *f);
PyObject *PConvIntVLAToPyList(int *f);

/* WARNING: the returned PyObject is unowned - it is intended for use
 * only for efficient detection of changes to dictionary values
 * following evaluation of some expression in the context of the
 * dictionary PAlter, PAlterState, etc. */
PyObject *PConvFloatToPyDictItem(PyObject *dict,char *key,float f); 
PyObject *PConvStringToPyDictItem(PyObject *dict,char *key,char *f);
PyObject *PConvIntToPyDictItem(PyObject *dict,char *key,int i);
/* end WARNING */

void PConvFloat3ToPyObjAttr(PyObject *obj,char *attr,float *v);
void PConvFloatToPyObjAttr(PyObject *obj,char *attr,float f);
void PConvIntToPyObjAttr(PyObject *obj,char *attr,int i);
void PConvInt2ToPyObjAttr(PyObject *obj,char *attr,int *v);
void PConvStringToPyObjAttr(PyObject *obj,char *attr,char *f);

int PConvPyObjectToFloat(PyObject *object,float *value);
int PConvPyObjectToInt(PyObject *object,int *value);

/* NOTE: the string routines will write strings up to the specified
 * length, PLUS a NULL...so watch out for array overruns */

int PConvPyObjectToStrMaxLen(PyObject *object,char *value,int ln);
int PConvPyObjectToStrMaxClean(PyObject *object,char *value,int ln);

PyObject *PConvStringListToPyList(int l,char **str);
PyObject *PConvStringVLAToPyList(char *str);

void PConv44PyListTo44f(PyObject *src,float *dest); /* note loss of precision */

int PConvPyListToFloatArray(PyObject *obj,float **f);
int PConvPyListToFloatArrayInPlace(PyObject *obj,float *ff,int ll);

PyObject *PConvFloatArrayToPyList(float *f,int l);

int PConvPyListToIntArray(PyObject *obj,int **f);
int PConvPyListToIntArrayInPlace(PyObject *obj,int *ff,int ll);

#endif







